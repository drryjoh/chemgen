#!/usr/bin/env python3

import os
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras import layers, callbacks, models, optimizers

# config
DATA_DIR        = "data_200"
MODEL_PATH      = "MLP_3.keras"
CSV_FILE        = "MLP_3.csv"
NUM_SAMPLES     = 200
M               = 96
FLAT_DIM        = M * M      
BATCH_SIZE      = 16
EPOCHS          = 100
HIDDEN_UNITS    = 16
LEARNING_RATE   = 1e-6       
CLIP_NORM       = 1.0        
VALIDATION_SPLIT= 0.2
RANDOM_SEED     = 42

# load and flatten data
X_list = []
# y_list = []  # comment out original inv_A loading
for i in range(NUM_SAMPLES):
    A  = np.loadtxt(os.path.join(DATA_DIR, f"A_{i}.csv"), delimiter=",", dtype=np.float32)
    # iA = np.loadtxt(os.path.join(DATA_DIR, f"inv_A_{i}.csv"), delimiter=",", dtype=np.float32)  # commented out
    X_list.append(A.ravel())
    # y_list.append(iA.ravel())

X = np.stack(X_list, axis=0)
# y = np.stack(y_list, axis=0)  # commented out original y

# load B and x data for new loss function
B_list, x_list = [], []
for i in range(NUM_SAMPLES):
    B = np.loadtxt(os.path.join(DATA_DIR, f"res_{i}.csv"), delimiter=",", dtype=np.float32)
    x = np.loadtxt(os.path.join(DATA_DIR, f"dy_{i}.csv"), delimiter=",", dtype=np.float32)
    B_list.append(B)
    x_list.append(x)
B = np.stack(B_list, axis=0)
x_true = np.stack(x_list, axis=0)

# prepare combined targets: [b, x]
y_combined = np.concatenate([B, x_true], axis=1)

# compute mean and std per feature for inputs only
eps    = 1e-7
X_mean = X.mean(axis=0)
X_std  = X.std(axis=0) + eps
# y_mean = y.mean(axis=0)  # commented out
# y_std  = y.std(axis=0) + eps
X_norm = (X - X_mean) / (X_std)
# y_norm = (y - y_mean) / (y_std)

# split into train and validation sets
X_tr, X_val, y_tr, y_val = train_test_split(
    X_norm, y_combined,
    test_size=VALIDATION_SPLIT,
    random_state=RANDOM_SEED,
    shuffle=True
)

###########
## MODEL ##
###########
model = tf.keras.Sequential([
    layers.Input(shape=(FLAT_DIM,)),
    layers.Dense(HIDDEN_UNITS, activation="relu"),
    layers.Dense(HIDDEN_UNITS, activation="relu"),
    layers.Dense(FLAT_DIM, activation="linear")
])

# compile model with gradient clipping
opt = optimizers.Adam(learning_rate=LEARNING_RATE, clipnorm=CLIP_NORM)
# model.compile(optimizer=opt, loss=tf.keras.losses.LogCosh())  # original loss commented out

# custom loss: ||P * b - x||^2
def custom_loss_fn(y_true, y_pred):
    # y_true: [b, x] concatenated
    b = y_true[:, :M]  # first M elements are b
    x = y_true[:, M:]  # next M elements are true x
    P = tf.reshape(y_pred, (-1, M, M))
    # predicted x = P @ b
    x_pred = tf.matmul(P, tf.expand_dims(b, -1))
    x_pred = tf.squeeze(x_pred, axis=-1)
    # mean squared error
    return tf.reduce_mean(tf.square(x_pred - x))

model.compile(optimizer=opt, loss=custom_loss_fn)

# set up early stopping and checkpoint
early_stop = callbacks.EarlyStopping(monitor="val_loss", patience=500, restore_best_weights=True)
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
print(f"final validation mse loss ||P*b - x||^2: {val_loss:.6f}")

# save normalization stats to csv for input only
with open(CSV_FILE, "w") as f:
    f.write("input_mean: ["  + ",".join(map(str, X_mean)) + "]\n")
    f.write("input_std:  ["  + ",".join(map(str, X_std )) + "]\n")
    # original output stats commented out
    # f.write("output_mean: [" + ",".join(map(str, y_mean)) + "]\n")
    # f.write("output_std:  [" + ",".join(map(str, y_std )) + "]\n")

model.summary()
