#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import scipy as sp
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras import layers, callbacks, models, optimizers

# config
DATA_DIR = "EULER_825"
MODEL_PATH = "bin/MLP_LU_1.keras"
CSV_FILE = "bin/MLP_LU_1.csv"
NUM_SAMPLES = 825
M = 96
FLAT_DIM = M * M
BATCH_SIZE = 16
EPOCHS = 10000
HIDDEN_UNITS = 4
LEARNING_RATE = 1e-5
CLIP_NORM = 0.1
VALIDATION_SPLIT = 0.3
RANDOM_SEED = 42
EPS = 1e-9

# load and flatten data
X_list, y_list = [], []
for i in range(NUM_SAMPLES):
    A = np.loadtxt(
        os.path.join(DATA_DIR, f"A_{i}.csv"), delimiter=",", dtype=np.float32
    )
    iA = np.loadtxt(
        os.path.join(DATA_DIR, f"A_inv_{i}.csv"), delimiter=",", dtype=np.float32
    )
    L, U = sp.linalg.lu(iA, permute_l=True)
    X_list.append(A.ravel())
    # y_list.append(iA.ravel())
    y_list.append(np.concatenate([L.ravel(), U.ravel()]))

X = np.stack(X_list, axis=0)
y = np.stack(y_list, axis=0)

# compute mean and std per feature

X_mean = X.mean(axis=0)
X_std = X.std(axis=0)
y_mean = y.mean(axis=0)
y_std = y.std(axis=0)

# normalize data
# X = (X - X_mean) / (X_std)
# y = (y - y_mean) / (y_std)

# split into train and validation sets
X_tr, X_val, y_tr, y_val = train_test_split(
    X, y, test_size=VALIDATION_SPLIT, random_state=RANDOM_SEED, shuffle=True
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
model = tf.keras.Sequential(
    [
        layers.Input(shape=(FLAT_DIM,)),
        layers.Rescaling(1.0 / (X_std + EPS), offset=-X_mean / (X_std + EPS)),

        # layers.Dense(HIDDEN_UNITS),
        # # layers.UnitNormalization(axis=-1),
        # layers.BatchNormalization(),
        # layers.Activation("gelu"),

        layers.Dense(HIDDEN_UNITS * 2),
        # layers.UnitNormalization(axis=-1),
        layers.BatchNormalization(),
        layers.Activation("gelu"),

        layers.Dense(HIDDEN_UNITS * 4),
        # layers.UnitNormalization(axis=-1),
        # layers.LayerNormalization(),
        layers.Activation("gelu"),

        layers.Dense(HIDDEN_UNITS * 8),
        # layers.UnitNormalization(axis=-1),
        # layers.LayerNormalization(),
        layers.Activation("gelu"),

        layers.Dense(FLAT_DIM * 2, activation="linear"),
        layers.Rescaling(y_std + EPS, offset=y_mean),
    ]
)
###########

# compile model with gradient clipping
opt = optimizers.Adam(learning_rate=LEARNING_RATE, clipnorm=CLIP_NORM)
model.compile(
    optimizer=opt,
    # loss=tf.keras.losses.MeanSquaredError()
    loss=tf.keras.losses.LogCosh()
    # loss=tf.keras.losses.MeanSquaredLogarithmicError()
    # loss=tf.keras.losses.MeanAbsolutePercentageError()
)

# set up early stopping and checkpoint
early_stop = callbacks.EarlyStopping(
    monitor="val_loss", patience=50, restore_best_weights=True
)
checkpoint = callbacks.ModelCheckpoint(MODEL_PATH, save_best_only=True)
# reduce_lr = callbacks.ReduceLROnPlateau(
#     monitor="val_loss", 
#     factor=0.25, 
#     patience=30 
#     min_lr=1e-5
# )

# train the model
history = model.fit(
    X_tr,
    y_tr,
    validation_data=(X_val, y_val),
    epochs=EPOCHS,
    batch_size=BATCH_SIZE,
    callbacks=[early_stop, 
               checkpoint, 
            #    reduce_lr
               ],
    verbose=1,
)

# evaluate on validation set
val_loss = model.evaluate(X_val, y_val, verbose=0)
print(f"final normalized-val mse: {val_loss:.6f}")

# # save normalization stats to csv
# with open(CSV_FILE, "w") as f:
#     f.write("input_mean: [" + ",".join(map(str, X_mean)) + "]\n")
#     f.write("input_std:  [" + ",".join(map(str, X_std)) + "]\n")
#     f.write("output_mean: [" + ",".join(map(str, y_mean)) + "]\n")
#     f.write("output_std:  [" + ",".join(map(str, y_std)) + "]\n")

model.summary()
