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
DATA_DIR        = "EULER_825"
FILE_NAME       = "CNN_1.keras"
MODEL_DIR       = "bin"
NUM_SAMPLES     = 825
M               = 96
FILTER_SIZE     = 8
KERNEL_SIZE     = 3
EPOCHS          = 50
BATCH_SIZE      = 16
LEARNING_RATE   = 1e-3
CLIP_NORM       = 0.1
VALIDATION_SPLIT= 0.3
RANDOM_SEED     = 42
EPS             = 1e-8

# Load and shape data
A_list = []
invA_list = []
for i in range(NUM_SAMPLES):
    A = np.loadtxt(os.path.join(DATA_DIR, f"A_{i}.csv"), delimiter=",", dtype=np.float32).reshape(M, M, 1)
    inv_A = np.loadtxt(os.path.join(DATA_DIR, f"A_inv_{i}.csv"), delimiter=",", dtype=np.float32).reshape(M, M)
    A_list.append(A)
    invA_list.append(inv_A)
A = np.stack(A_list, axis=0)         
inv_A = np.stack(invA_list, axis=0)  

# Normalize data
A_mean, A_std = A.mean(axis=0), A.std(axis=0) + EPS
A_inv_std, A_inv_mean = inv_A.std(axis=0) + EPS, inv_A.mean(axis=0)

# Shuffle and split data
dataset = tf.data.Dataset.from_tensor_slices((A, inv_A)).shuffle(512, seed=RANDOM_SEED, reshuffle_each_iteration=True)
val_size = int(VALIDATION_SPLIT * NUM_SAMPLES)
train_ds = dataset.skip(val_size).batch(BATCH_SIZE).prefetch(tf.data.AUTOTUNE)
val_ds   = dataset.take(val_size).batch(BATCH_SIZE).prefetch(tf.data.AUTOTUNE)

# Define CNN model
model = models.Sequential([
    layers.Input(shape=(M, M, 1)),
    layers.Rescaling(1.0/A_std, offset=-A_mean/A_std),

    layers.Conv2D(FILTER_SIZE, KERNEL_SIZE, padding='same', activation='gelu'),
    # layers.LayerNormalization(axis=[1,2,3]),
    layers.Conv2D(FILTER_SIZE, KERNEL_SIZE, padding='same', activation='gelu'),
    layers.Conv2D(FILTER_SIZE, KERNEL_SIZE, padding='same', activation='gelu'),
    layers.Conv2D(1, KERNEL_SIZE, padding='same'),
    layers.Reshape((96,96)),

    layers.Rescaling(A_inv_std, offset=A_inv_mean)
])

# Compile with MSE loss
opt = optimizers.Adam(learning_rate=LEARNING_RATE, clipnorm=CLIP_NORM)
model.compile(optimizer=opt,
            #   loss="log_cosh", 
              loss='mse'
            #   loss=tf.keras.losses.MeanSquaredLogarithmicError()
            # loss=tf.keras.losses.MeanAbsolutePercentageError()
              )

# Callbacks
early_stop = callbacks.EarlyStopping(monitor="val_loss", patience=30, restore_best_weights=True)
checkpoint = callbacks.ModelCheckpoint(os.path.join(MODEL_DIR, FILE_NAME), save_best_only=True)

# Train
history = model.fit(
    train_ds,
    validation_data=val_ds,
    epochs=EPOCHS,
    callbacks=[early_stop, checkpoint],
    verbose=2
)

# Evaluate
val_loss = model.evaluate(val_ds, verbose=0)
print(f"Final validation MSE: {val_loss:.6f}")

# Model Summary
model.summary()
