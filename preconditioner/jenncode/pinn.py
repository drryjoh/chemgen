import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# ------------------------------------------------------------------
# 1. CNN architecture that takes an n×n matrix (here 96×96) and
#    returns another n×n matrix: an approximation of A^{-1}.
# ------------------------------------------------------------------
model = keras.Sequential([
    layers.Conv2D(16, 3, padding='same', activation='relu', input_shape=(96, 96, 1)),
    layers.Conv2D( 8, 3, padding='same', activation='relu'),
    layers.Conv2D( 1, 3, padding='same', activation=None),   # linear output
    layers.Reshape((96, 96, 1))                              # make sure shape matches
])

optimizer   = tf.keras.optimizers.Adam(learning_rate=1e-4)
batch_size  = 1
epochs      = 20_000

# ------------------------------------------------------------------
# 2. Build a tf.data pipeline that feeds only (A, b) pairs.
#    Shapes:  A : (N, 96, 96, 1)   b : (N, 96, 1)
# ------------------------------------------------------------------
train_ds = tf.data.Dataset.from_tensor_slices((A, b)) \
                           .shuffle(buffer_size=len(A)) \
                           .batch(batch_size, drop_remainder=True)

# ------------------------------------------------------------------
# 3. Training loop: residual-only loss  ||A xθ − b||²
# ------------------------------------------------------------------
for epoch in range(epochs):
    epoch_loss = 0.0

    for step, (A_batch, b_batch) in enumerate(train_ds):
        with tf.GradientTape() as tape:

            # ----- forward pass -----------------------------------
            A_inv_pred = model(A_batch, training=True)               # (B,96,96,1)
            A_inv_mat  = tf.reshape(A_inv_pred, [-1, 96, 96])        # (B,96,96)
            b_vec      = tf.reshape(b_batch,   [-1, 96, 1])          # (B,96,1)
            x_pred     = tf.matmul(A_inv_mat, b_vec)                 # (B,96,1)

            # ----- residual  r = Axθ − b --------------------------
            A_mat      = tf.reshape(A_batch,  [-1, 96, 96])
            residual   = tf.matmul(A_mat, x_pred) - b_vec            # (B,96,1)
            loss       = tf.reduce_mean(tf.square(residual))         # scalar

        # ----- back-prop & update ---------------------------------
        grads = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(grads, model.trainable_variables))
        epoch_loss += loss.numpy()

    # average minibatch loss for the epoch
    epoch_loss /= (step + 1)

    if (epoch + 1) % 100 == 0:
        print(f'Epoch {epoch+1:05d} | residual MSE {epoch_loss:.6e}")

# ------------------------------------------------------------------
# 4. Save the trained inverse operator
# ------------------------------------------------------------------
model.save("cnn_pinn_Ainv.keras")
