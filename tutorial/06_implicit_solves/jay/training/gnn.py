script = '''
"""
Example Sequential GNN with Spektral + TensorFlow Keras
------------------------------------------------------

• GraphConv  (Kipf & Welling GCN)
• GraphSageConv (Hamilton et al. GraphSAGE)
• GATConv   (Veličković et al. Graph Attention Network)
• GlobalSumPool  (readout)
------------------------------------------------------

This minimal script builds the model, creates random
graph data, compiles, and trains for a single epoch.
Feel free to copy‑paste into a notebook or .py file.
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, Sequential
from spektral.layers import (
    GraphConv,
    GraphSageConv,
    GATConv,
    GlobalSumPool,
)

# ---- Random toy dataset -------------------------------------------------
N         = 100          # number of nodes
F         = 16           # node feature dimensionality
n_classes = 3            # for soft‑max classification

# Node features
X = np.random.randn(N, F).astype("float32")

# Adjacency matrix (undirected, no self‑loops)
A = np.random.randint(0, 2, (N, N)).astype("float32")
A = np.maximum(A, A.T)          # make symmetric
np.fill_diagonal(A, 0)          # remove self‑loops

# Graph‑level label
y = tf.keras.utils.to_categorical(
        np.random.randint(0, n_classes, 1), num_classes=n_classes
    ).astype("float32")

# ---- Sequential GNN model ----------------------------------------------
model = Sequential(
    [
        GraphConv(32, activation="relu"),
        GraphSageConv(32, activation="relu"),
        GATConv(32, activation="elu", heads=4, concat=False),
        GlobalSumPool(),
        layers.Dense(64, activation="relu"),
        layers.Dense(n_classes, activation="softmax"),
    ],
    name="GNN_sequential_demo",
)

model.compile(optimizer="adam", loss="categorical_crossentropy", metrics=["acc"])
model.summary()

# ---- Dummy train --------------------------------------------------------
# Spektral layers expect inputs as ([X, A], y) tuples for single graphs.
model.fit(x=[X, A], y=y, epochs=1, verbose=2)
'''

print(script)
