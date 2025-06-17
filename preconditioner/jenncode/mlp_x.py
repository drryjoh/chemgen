#!/usr/bin/env python3

import os
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras import layers, callbacks, optimizers

# config
DATA_DIR         = "EULER_825"
MODEL_PATH       = "MLP_x_1.keras"
CSV_FILE         = "MLP_x_1.csv"
NUM_SAMPLES      = 825
M                = 96
FLAT_DIM         = M * M      
BATCH_SIZE       = 64
EPOCHS           = 20000
HIDDEN_UNITS     = 8
LEARNING_RATE    = 1e-6
CLIP_NORM        = 1.0        
VALIDATION_SPLIT = 0.3
RANDOM_SEED      = 42
EPS              = 1e-8

# 1) load and flatten input matrices a into X_list
X_list = []
for i in range(NUM_SAMPLES):
    A = np.loadtxt(os.path.join(DATA_DIR, f"A_{i}.csv"),
                   delimiter=",", dtype=np.float32)
    X_list.append(A.ravel())
X = np.stack(X_list, axis=0)

# 2) load b vectors and true x vectors
B_list, x_list = [], []
for i in range(NUM_SAMPLES):
    B = np.loadtxt(os.path.join(DATA_DIR, f"res_{i}.csv"),
                   delimiter=",", dtype=np.float32)
    x = np.loadtxt(os.path.join(DATA_DIR, f"dy_{i}.csv"),
                   delimiter=",", dtype=np.float32)
    B_list.append(B)
    x_list.append(x)
B = np.stack(B_list, axis=0)       # shape: (NUM_SAMPLES, M)
x_true = np.stack(x_list, axis=0)  # shape: (NUM_SAMPLES, M)

# 3) compute normalization statistics
X_mean = X.mean(axis=0)
X_std  = X.std(axis=0) + EPS
B_mean = B.mean(axis=0)
B_std  = B.std(axis=0) + EPS
x_mean = x_true.mean(axis=0)
x_std  = x_true.std(axis=0) + EPS

# 4) normalize inputs and targets
X_norm = (X - X_mean) / X_std
B_norm = (B - B_mean) / B_std
x_norm = (x_true - x_mean) / x_std

# 5) prepare combined normalized targets [b_norm, x_norm]
y_norm = np.concatenate([B_norm, x_norm], axis=1)

# 6) split into train and validation sets
X_tr, X_val, y_tr, y_val = train_test_split(
    X_norm, y_norm,
    test_size=VALIDATION_SPLIT,
    random_state=RANDOM_SEED,
    shuffle=True
)

# 7) build the model
model = tf.keras.Sequential([
    layers.Input(shape=(FLAT_DIM,)),
    layers.Dense(HIDDEN_UNITS),
    layers.BatchNormalization(),
    layers.LeakyReLU(negative_slope=0.1),
    layers.Dense(HIDDEN_UNITS),
    layers.LeakyReLU(negative_slope=0.1),
    layers.Dense(HIDDEN_UNITS),
    layers.LeakyReLU(negative_slope=0.1),
    layers.Dense(FLAT_DIM, activation="linear")  # predicts flattened P
])

# custom loss: mse between normalized x_pred and x_true_norm
def custom_loss_fn(y_true, y_pred):
    b_norm       = y_true[:, :M]      # first M entries are normalized B
    x_true_norm  = y_true[:, M:]      # next M entries are normalized x
    P            = tf.reshape(y_pred, (-1, M, M))
    x_pred_norm  = tf.matmul(P, tf.expand_dims(b_norm, -1))
    x_pred_norm  = tf.squeeze(x_pred_norm, axis=-1)
    return tf.reduce_mean(tf.square(x_pred_norm - x_true_norm))

# compile with adam optimizer and gradient clipping
opt = optimizers.Adam(learning_rate=LEARNING_RATE, clipnorm=CLIP_NORM)
model.compile(optimizer=opt, loss=custom_loss_fn)

# callbacks for early stopping and checkpointing
early_stop = callbacks.EarlyStopping(
    monitor="val_loss", patience=500, restore_best_weights=True
)
checkpoint = callbacks.ModelCheckpoint(MODEL_PATH, save_best_only=True)

# train the model
history = model.fit(
    X_tr, y_tr,
    validation_data=(X_val, y_val),
    epochs=EPOCHS,
    batch_size=BATCH_SIZE,
    callbacks=[early_stop, checkpoint],
    verbose=1
)

# evaluate on validation set
val_loss = model.evaluate(X_val, y_val, verbose=0)
print(f"\n\nfinal validation mse loss ||p*b - x||^2: {val_loss:.6f}")

# compute output (p) normalization stats and save only input/output stats
y_pred_norm = model.predict(X_tr, batch_size=BATCH_SIZE)
output_mean = y_pred_norm.mean(axis=0)
output_std  = y_pred_norm.std(axis=0) + EPS

with open(CSV_FILE, "w") as f:
    # save input normalization
    f.write("input_mean:  [" + ",".join(map(str, X_mean))    + "]\n")
    f.write("input_std:   [" + ",".join(map(str, X_std))     + "]\n")
    # save output normalization
    f.write("output_mean: [" + ",".join(map(str, output_mean)) + "]\n")
    f.write("output_std:  [" + ",".join(map(str, output_std))  + "]\n")

# # example unnormalize p for a new matrix a_new
# a_flat       = a_new.ravel()
# a_norm       = (a_flat - X_mean) / X_std
# p_norm_flat  = model.predict(a_norm[None, :])[0]
# p_flat       = p_norm_flat * output_std + output_mean
# p            = p_flat.reshape(M, M)

# model summary
model.summary()
