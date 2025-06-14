# import numpy as np
# import tensorflow as tf
# from tensorflow import keras
# from tensorflow.keras import layers, models, optimizers, losses, callbacks
# # np.set_printoptions(threshold=np.inf)

# # load data
# A_inv = np.loadtxt('./data_1/A_inv.csv', delimiter=',').astype(np.float32)
# A = np.loadtxt('./data_1/A.csv', delimiter=',').astype(np.float32)
# x = np.loadtxt('./data_1/dy.csv', delimiter=',').astype(np.float32)
# b = np.loadtxt('./data_1/res.csv', delimiter=',').astype(np.float32)
# # print(A_insv)

# # rearrange shape
# # A_inv_b = A_inv.reshape(-1, 96, 96, 1)
# A_inv = A_inv.reshape(-1,96,96,1)
# A = A.reshape(-1, 96, 96, 1)
# x = x.reshape(-1, 96, 1)
# b = b.reshape(-1, 96, 1)

# # calculate inverse A * b
# # A_inv_b = np.linalg.matmul(A_inv, b)
# # print(A_inv_b)

# # build sequential cnn
# model = keras.Sequential([
#     layers.Conv2D(16, 3, padding='same', activation='relu', input_shape=(96, 96, 1)),
#     # layers.Conv2D(16, 3, padding='same', activation='relu'),
#     layers.Conv2D( 8, 3, padding='same', activation='relu'),
#     layers.Conv2D( 1, 3, padding='same', activation=None),
#     layers.Reshape((96, 96, 1))  
# ])
# # model = keras.Sequential([
# #     layers.SeparableConv2D(16, 3, padding='same', activation='relu', input_shape=(96, 96, 1)),
# #     layers.SeparableConv2D(8, 3, padding='same', activation='relu'),
# #     layers.SeparableConv2D(1, 3, padding='same', activation=None),
# #     layers.Reshape((96, 96))
# # ])

# optimizer = tf.keras.optimizers.Adam(learning_rate=1e-4)

# # training loop: iterate through epochs and batches
# batch_size = 1
# epochs = 20000
# how_many_datasets = A.shape[0]
# steps_per_epoch = how_many_datasets // batch_size

# for epoch in range(epochs):
#     # randomly shuffle datasets (indices) for each epoch 
#     # idx = np.random.permutation(how_many_datasets)
#     # A_shuf = A[idx]
#     # b_shuf = b[idx]
#     # x_shuf = x[idx]
#     epoch_loss = 0.0

#     for step in range(steps_per_epoch):
#         # calculate batch start and end indices
#         start = step * batch_size
#         end = start + batch_size
#         # convert batch data to tensors
#         A_batch = tf.convert_to_tensor(A[start:end])
#         b_batch = tf.convert_to_tensor(b[start:end])
#         x_batch = tf.convert_to_tensor(x[start:end])
#         A_inv_batch = tf.convert_to_tensor(A_inv[start:end])
#         # A_inv_b_batch = tf.convert_to_tensor(A_inv_b[start:end])

#         with tf.GradientTape() as tape:

#             # forward pass: predict p for the current batch
#             # P_pred = model(A_batch, training=True)
#             A_inv_pred =  model(A_batch, training=True)

#             # compute predicted x by multiplying p and b
#             # P_b = tf.linalg.solve(P_pred, b_batch)
#             A_inv_flat = tf.reshape(A_inv_pred, [-1, 96, 96])      # (1,96,96)
#             b_flat      = tf.reshape(b_batch,    [-1, 96, 1])      # (1,96,1)
#             x_pred      = tf.matmul(A_inv_flat,  b_flat)  

#             # calculate mean squared error loss between predicted and true x
#             loss = tf.reduce_mean(tf.square(x_pred - x_batch))

#         # compute gradients and update weights
#         grads = tape.gradient(loss, model.trainable_variables)
#         optimizer.apply_gradients(zip(grads, model.trainable_variables))
#         epoch_loss += loss.numpy()

#     # average loss for the epoch
#     epoch_loss /= steps_per_epoch
#     # print epoch number and corresponding loss
#     print(f'epoch {epoch+1:03d}, loss {epoch_loss:.6f}')

# # save the trained model
# model.save('cnn_3.keras')









































































#!/usr/bin/env python3

import os
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, callbacks, models, optimizers

# Config
DATA_DIR        = "data_200"
MODEL_PATH      = "CNN_1.keras"
NUM_SAMPLES     = 200
M               = 96
EPOCHS          = 10000
BATCH_SIZE      = 32
LEARNING_RATE   = 1e-6
CLIP_NORM       = 1.0
VALIDATION_SPLIT= 0.2
RANDOM_SEED     = 42

# Load and shape data
A_list = []
invA_list = []
for i in range(NUM_SAMPLES):
    A = np.loadtxt(os.path.join(DATA_DIR, f"A_{i}.csv"), delimiter=",", dtype=np.float32).reshape(M, M, 1)
    inv_A = np.loadtxt(os.path.join(DATA_DIR, f"inv_A_{i}.csv"), delimiter=",", dtype=np.float32).reshape(M, M)
    A_list.append(A)
    invA_list.append(inv_A)

A = np.stack(A_list, axis=0)          # Shape: (200, 96, 96, 1)
inv_A = np.stack(invA_list, axis=0)   # Shape: (200, 96, 96)

# Shuffle and split data
dataset = tf.data.Dataset.from_tensor_slices((A, inv_A)).shuffle(512, seed=RANDOM_SEED, reshuffle_each_iteration=True)
val_size = int(VALIDATION_SPLIT * NUM_SAMPLES)
train_ds = dataset.skip(val_size).batch(BATCH_SIZE).prefetch(tf.data.AUTOTUNE)
val_ds   = dataset.take(val_size).batch(BATCH_SIZE).prefetch(tf.data.AUTOTUNE)

# Define CNN model
model = models.Sequential([
    layers.Input(shape=(M, M, 1)),
    layers.Conv2D(16, kernel_size=3, padding='same', activation='relu'),
    layers.Conv2D(8, kernel_size=3, padding='same', activation='relu'),
    layers.Conv2D(1, kernel_size=3, padding='same'),
    layers.Reshape((M, M))
])

# Compile with MSE loss
opt = optimizers.Adam(learning_rate=LEARNING_RATE, clipnorm=CLIP_NORM)
model.compile(optimizer=opt, loss='mse')

# Callbacks
early_stop = callbacks.EarlyStopping(monitor="val_loss", patience=30, restore_best_weights=True)
checkpoint = callbacks.ModelCheckpoint(MODEL_PATH, save_best_only=True)

# Train
history = model.fit(
    train_ds,
    validation_data=val_ds,
    epochs=EPOCHS,
    callbacks=[early_stop, checkpoint],
    verbose=1
)

# Evaluate
val_loss = model.evaluate(val_ds, verbose=0)
print(f"Final validation MSE: {val_loss:.6f}")

# Model Summary
model.summary()
