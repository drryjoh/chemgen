#!/usr/bin/env python3

import os
import numpy as np
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras import layers, callbacks, optimizers
from tensorflow.keras.optimizers import Adam


# config
DATA_DIR         = "ROSENBROCK_400"
MODEL_PATH       = "bin/MLP_x_1_R.keras"
CSV_FILE         = "bin/MLP_x_1_R.csv"
NUM_SAMPLES      = 192
M                = 96
FLAT_DIM         = M * M      
BATCH_SIZE       = 4
EPOCHS           = 20000
HIDDEN_UNITS     = 8
LEARNING_RATE    = 1e-6
CLIP_NORM        = 1.0        
VALIDATION_SPLIT = 0.5
RANDOM_SEED      = 42
EPS              = 1e-12

#load A and A_inv
X_list, A_inv_list = [], []
for i in range(NUM_SAMPLES):
    A = np.loadtxt(os.path.join(DATA_DIR, f"A_{i}.csv"),
                   delimiter=",", dtype=np.float32)
    X_list.append(A.ravel())
    A_inv = np.loadtxt(os.path.join(DATA_DIR, f"A_inv_{i}.csv"),
                   delimiter=",", dtype=np.float32)
    A_inv_list.append(A_inv.ravel())
A = np.stack(X_list, axis=0)
A_inv = np.stack(A_inv_list, axis=0)

#load b and x
b_list, x_list = [], []
for i in range(NUM_SAMPLES):
    b = np.loadtxt(os.path.join(DATA_DIR, f"b_{i}.csv"),
                   delimiter=",", dtype=np.float32)
    x = np.loadtxt(os.path.join(DATA_DIR, f"x_{i}.csv"),
                   delimiter=",", dtype=np.float32)
    b_list.append(b)
    x_list.append(x)
b = np.stack(b_list, axis=0)      
x = np.stack(x_list, axis=0)  

#normalize A, b, x, and A_inv
A_mean = A.mean(axis=0)
A_std  = A.std(axis=0) + EPS
b_mean = b.mean(axis=0)
b_std  = b.std(axis=0) + EPS
x_mean = x.mean(axis=0)
x_std  = x.std(axis=0) + EPS
A_inv_mean = A_inv.mean(axis=0)
A_inv_std = A_inv.std(axis=0) + EPS
###
A_norm = (A - A_mean) / A_std
b_norm = (b - b_mean) / b_std
x_norm = (x - x_mean) / x_std
A_inv_norm = (A_inv - A_inv_mean) / A_inv_std

#merge target training target
# y_norm = np.concatenate([b_norm, x_norm], axis=1)

#split to training and validation set
A_tr, A_val, b_tr, b_val, A_inv_tr, A_inv_val = train_test_split(
    # X_norm, y_norm,
    A_norm, b_norm, A_inv_norm,
    test_size=VALIDATION_SPLIT,
    random_state=RANDOM_SEED,
    shuffle=True
)

#convert to data sets
train_ds = tf.data.Dataset.from_tensor_slices((A_tr, b_tr, A_inv_tr))
train_ds = train_ds.shuffle(buffer_size=1024).batch(BATCH_SIZE)
val_ds = tf.data.Dataset.from_tensor_slices((A_val, b_val, A_inv_val))
val_ds = val_ds.batch(BATCH_SIZE)

# build model
model = tf.keras.Sequential([
    layers.Input(shape=(FLAT_DIM,)),
    layers.Dense(HIDDEN_UNITS), layers.BatchNormalization(), layers.LeakyReLU(negative_slope=0.001),
    layers.Dense(HIDDEN_UNITS), layers.LeakyReLU(negative_slope=0.001),
    layers.Dense(HIDDEN_UNITS), layers.LeakyReLU(negative_slope=0.001),
    layers.Dense(FLAT_DIM, activation="linear")
])

#######################
## for ||P*b - x||^2 ##
#######################
# def custom_loss_fn(y_true, y_pred):
#     b_norm       = y_true[:, :M]     
#     x_true_norm  = y_true[:, M:]     
#     P            = tf.reshape(y_pred, (-1, M, M))
#     x_pred_norm  = tf.matmul(P, tf.expand_dims(b_norm, -1))
#     x_pred_norm  = tf.squeeze(x_pred_norm, axis=-1)
#     return tf.reduce_mean(tf.square(x_pred_norm - x_true_norm))

# compile with adam optimizer and gradient clipping
# opt = optimizers.Adam(learning_rate=LEARNING_RATE, clipnorm=CLIP_NORM)
# model.compile(optimizer=opt, loss=custom_loss_fn)

# # callbacks for early stopping and checkpointing
# early_stop = callbacks.EarlyStopping(
#     monitor="val_loss", patience=500, restore_best_weights=True
# )
# checkpoint = callbacks.ModelCheckpoint(MODEL_PATH, save_best_only=True)
# reduce_lr = callbacks.ReduceLROnPlateau(
#     monitor="val_loss", 
#     factor=0.025, 
#     patience=100, 
#     min_lr=1e-10
# )

# # train the model
# history = model.fit(
#     X_tr, y_tr,
#     validation_data=(X_val, y_val),
#     epochs=EPOCHS,
#     batch_size=BATCH_SIZE,
#     callbacks=[early_stop, checkpoint, reduce_lr],
#     verbose=2
# )

# # evaluate on validation set
# val_loss = model.evaluate(X_val, y_val, verbose=0)
# print(f"\n\nfinal validation mse loss ||p*b - x||^2: {val_loss:.6f}")

# # compute output (p) normalization stats and save only input/output stats
# y_pred_norm = model.predict(X_tr, batch_size=BATCH_SIZE)
# output_mean = y_pred_norm.mean(axis=0)
# output_std  = y_pred_norm.std(axis=0) + EPS
###
##
#

#############################
## for ||P*b - A_inv*b||^2 ##
#############################
opt = optimizers.Adam(learning_rate=LEARNING_RATE, clipnorm=CLIP_NORM)

def custom_loss_fn(P_flat, b, Ainv):
    P = tf.reshape(P_flat, (-1, M, M)) 
    Ainv = tf.reshape(Ainv, (-1, M, M))         
    b = tf.expand_dims(b, -1)               
    Pb = tf.matmul(P, b)                   
    Ainv_b = tf.matmul(Ainv, b)            
    Pb = tf.squeeze(Pb, axis=-1)
    Ainv_b = tf.squeeze(Ainv_b, axis=-1)
    return tf.reduce_mean(tf.square(Pb - Ainv_b))

@tf.function
def train_step(x, b, Ainv):
    with tf.GradientTape() as tape:
        P_flat = model(x, training=True)
        loss = custom_loss_fn(P_flat, b, Ainv)
    grads = tape.gradient(loss, model.trainable_variables)
    opt.apply_gradients(zip(grads, model.trainable_variables))
    return loss

@tf.function
def val_step(x, b, Ainv):
    P_flat = model(x, training=False)
    return custom_loss_fn(P_flat, b, Ainv)

#custom training
best_val_loss = float('inf')
patience_counter = 0
for epoch in range(EPOCHS):

    #training
    train_loss = 0.0
    for x_batch, b_batch, Ainv_batch in train_ds:
        train_loss += train_step(x_batch, b_batch, Ainv_batch).numpy()
    train_loss /= len(train_ds)

    #validation
    val_loss = 0.0
    for x_batch, b_batch, Ainv_batch in val_ds:
        val_loss += val_step(x_batch, b_batch, Ainv_batch).numpy()
    val_loss /= len(val_ds)

    print(f"Epoch {epoch+1:4d} | Train Loss: {train_loss:.6f} | Val Loss: {val_loss:.6f}")

    #early stopping
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        patience_counter = 0
        model.save(MODEL_PATH)
    else:
        patience_counter += 1
        if patience_counter >= 500:
            print("Early stopping triggered.")
            break

all_preds = []
for x_batch, _, _ in train_ds:
    pred = model(x_batch, training=False)
    all_preds.append(pred.numpy())

P_preds = np.concatenate(all_preds, axis=0)  
output_mean = P_preds.mean(axis=0)
output_std = P_preds.std(axis=0) + EPS
###
##
#

with open(CSV_FILE, "w") as f:
    f.write("input_mean:  [" + ",".join(map(str, A_mean)) + "]\n")
    f.write("input_std:   [" + ",".join(map(str, A_std)) + "]\n")
    f.write("output_mean: [" + ",".join(map(str, output_mean)) + "]\n")
    f.write("output_std:  [" + ",".join(map(str, output_std)) + "]\n")

# # example unnormalize p for a new matrix a_new
# a_flat       = a_new.ravel()
# a_norm       = (a_flat - X_mean) / X_std
# p_norm_flat  = model.predict(a_norm[None, :])[0]
# p_flat       = p_norm_flat * output_std + output_mean
# p            = p_flat.reshape(M, M)

# model summary
model.summary()
