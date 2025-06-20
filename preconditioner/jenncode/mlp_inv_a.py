# import numpy as np
# import tensorflow as tf

# A = np.loadtxt('./data_1/A.csv', delimiter=',').astype('float32').reshape(-1,96,96,1)
# b = np.loadtxt('./data_1/res.csv', delimiter=',').astype('float32').reshape(-1,96,1)
# x = np.loadtxt('./data_1/dy.csv', delimiter=',').astype('float32').reshape(-1,96,1)
# y = np.concatenate([b, x], axis=-1)

# ds = tf.data.Dataset.from_tensor_slices((A, y)).shuffle(512, reshuffle_each_iteration=True)
# val_size = int(0.1 * A.shape[0])
# train_ds = ds.skip(val_size).batch(32).prefetch(tf.data.AUTOTUNE)
# val_ds   = ds.take(val_size).batch(32).prefetch(tf.data.AUTOTUNE)

# model = tf.keras.Sequential([
#   tf.keras.Input(shape=(96,96,1)),
#   tf.keras.layers.Dense(8, activation='relu'),
#   tf.keras.layers.Dense(4, activation='relu'),
#   tf.keras.layers.Dense(1, activation='linear'),     # ← 1 channel per pixel
#   tf.keras.layers.Reshape((96, 96))                  # now reshapes 96×96×1 → 96×96
# ])


# def custom_loss(y_true, y_pred):
#     P     = y_pred
#     b_vec = y_true[..., 0]
#     x_vec = y_true[..., 1]
#     y_hat = tf.linalg.matvec(P, b_vec)
#     return tf.reduce_mean(tf.square(y_hat - x_vec))

# model.compile(
#     optimizer=tf.keras.optimizers.Adam(1e-4),
#     loss=custom_loss
# )

# es  = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=30, restore_best_weights=True)
# rlr = tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=10, min_lr=1e-7)

# model.fit(train_ds, validation_data=val_ds, epochs=11000, callbacks=[es, rlr], verbose=2)
# model.save('mlp_1.keras')


































































#!/usr/bin/env python3

import os
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras import layers, callbacks, models, optimizers

# config
DATA_DIR        = "EULER_825"
MODEL_PATH      = "TEST.keras"
CSV_FILE        = "TEST.csv"
NUM_SAMPLES     = 825
M               = 96
FLAT_DIM        = M * M      
BATCH_SIZE      = 16
EPOCHS          = 100
HIDDEN_UNITS    = 8
LEARNING_RATE   = 1e-3
CLIP_NORM       = 1.0        
VALIDATION_SPLIT= 0.3
RANDOM_SEED     = 42

# load and flatten data
X_list, y_list = [], []
for i in range(NUM_SAMPLES):
    A  = np.loadtxt(os.path.join(DATA_DIR, f"A_{i}.csv"), delimiter=",", dtype=np.float32)
    iA = np.loadtxt(os.path.join(DATA_DIR, f"A_inv_{i}.csv"), delimiter=",", dtype=np.float32)
    X_list.append(A.ravel())
    y_list.append(iA.ravel())

X = np.stack(X_list, axis=0)
y = np.stack(y_list, axis=0)

# compute mean and std per feature
eps    = 1e-8
X_mean = X.mean(axis=0)
X_std  = X.std(axis=0) + eps
y_mean = y.mean(axis=0)
y_std  = y.std(axis=0) + eps
X_norm = (X - X_mean) / (X_std)
y_norm = (y - y_mean) / (y_std)

# split into train and validation sets
X_tr, X_val, y_tr, y_val = train_test_split(
    X_norm, y_norm,
    test_size=VALIDATION_SPLIT,
    random_state=RANDOM_SEED,
    shuffle=True
)

###########
## MODEL ##
###########
# inputs = layers.Input(shape=(FLAT_DIM,))
# x = layers.Dense(HIDDEN_UNITS, activation=None)(inputs)
# x = layers.BatchNormalization()(x) 
# x = layers.Activation("relu")(x)
# x = layers.Dense(HIDDEN_UNITS, activation=None)(x)
# x = layers.BatchNormalization()(x)
# x = layers.Activation("relu")(x)
# x = layers.Dense(HIDDEN_UNITS, activation=None)(x)
# x = layers.BatchNormalization()(x)
# x = layers.Activation("relu")(x)
# outputs = layers.Dense(FLAT_DIM, activation=None)(x)
# model = models.Model(inputs, outputs)
###########
model = tf.keras.Sequential([
    layers.Input(shape=(FLAT_DIM,)),
    layers.Dense(HIDDEN_UNITS),
    # layers.UnitNormalization(axis=-1),
    layers.BatchNormalization(),
    # layers.Activation("relu"),
    layers.LeakyReLU(negative_slope=0.01),
    layers.Dense(HIDDEN_UNITS),
    # layers.BatchNormalization(),
    # layers.Activation("relu"),
    layers.LeakyReLU(negative_slope=0.01),
    layers.Dense(FLAT_DIM, activation="linear")
])
###########

# compile model with gradient clipping
opt = optimizers.Adam(learning_rate=LEARNING_RATE, clipnorm=CLIP_NORM)
model.compile(optimizer=opt, 
            #   loss=tf.keras.losses.MeanSquaredError()
              loss=tf.keras.losses.LogCosh()
            #   loss='binary_crossentropy',
            #   metrics=['accuracy']
              )

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
print(f"final normalized-val mse: {val_loss:.6f}")

# save normalization stats to csv
with open(CSV_FILE, "w") as f:
    f.write("input_mean: ["  + ",".join(map(str, X_mean)) + "]\n")
    f.write("input_std:  ["  + ",".join(map(str, X_std )) + "]\n")
    f.write("output_mean: [" + ",".join(map(str, y_mean)) + "]\n")
    f.write("output_std:  [" + ",".join(map(str, y_std )) + "]\n")

model.summary()