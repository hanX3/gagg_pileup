"""Versioned model-family routing; no training data is opened for prediction."""
from .alpha_pair_ml import blind_predict
from .alpha_pair_study import load_model as load_unordered,PAIR_VARIANTS
from .alpha_pair_ordered import load_model as load_ordered,predict as ordered_predict,VARIANTS as ORDERED

ALL_PAIRS=PAIR_VARIANTS+ORDERED


def load_model(name):
    return load_ordered(name) if name.startswith('pair_ordered_') else load_unordered(name)


def predict(model,voltage,threshold):
    return ordered_predict(model,voltage,threshold) if hasattr(model,'count_head') else blind_predict(model,voltage,threshold=threshold)
