#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict


def read_features(path):
    df = pd.read_csv(path, sep="\t")
    if "Sample" not in df.columns:
        raise ValueError(f"{path} must contain a Sample column")
    return df.set_index("Sample")


def read_metadata(path):
    meta = pd.read_csv(path, sep="\t")
    meta["Sample"] = meta["Sample"].astype(str).str.replace(".t2", "-t2", regex=False)
    meta = meta[meta["Group"].isin(["HC", "SCZ"])].copy()
    return meta.set_index("Sample")


def run_model(feature_path, metadata_path, output_prefix, repeats=10, folds=10, seed=2024):
    X = read_features(feature_path)
    meta = read_metadata(metadata_path)
    common = X.index.intersection(meta.index)
    X = X.loc[common].fillna(0)
    y = meta.loc[common, "Group"].astype(str)
    y_binary = (y == "SCZ").astype(int)

    cv = RepeatedStratifiedKFold(n_splits=folds, n_repeats=repeats, random_state=seed)
    model = RandomForestClassifier(n_estimators=1000, random_state=seed, n_jobs=-1)
    pred = cross_val_predict(model, X, y_binary, cv=cv, method="predict_proba", n_jobs=-1)[:, 1]
    auc = roc_auc_score(y_binary, pred)

    model.fit(X, y_binary)
    importance = pd.DataFrame({
        "feature": X.columns,
        "importance": model.feature_importances_
    }).sort_values("importance", ascending=False)

    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"Sample": common, "Group": y.values, "SCZ_probability": pred}).to_csv(
        str(output_prefix) + ".crossval_predictions.tsv", sep="\t", index=False
    )
    importance.to_csv(str(output_prefix) + ".feature_importance.tsv", sep="\t", index=False)
    pd.DataFrame({"model": [output_prefix.name], "auc": [auc], "n_samples": [len(common)], "n_features": [X.shape[1]]}).to_csv(
        str(output_prefix) + ".auc.tsv", sep="\t", index=False
    )


def main():
    parser = argparse.ArgumentParser(description="Random-forest analysis for SCZ vs HC microbial feature matrices.")
    parser.add_argument("-i", "--input", required=True, help="Feature matrix with samples in rows and a Sample column.")
    parser.add_argument("-m", "--metadata", required=True, help="Sample metadata with Sample and Group columns.")
    parser.add_argument("-o", "--output-prefix", required=True, help="Output prefix.")
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--folds", type=int, default=10)
    parser.add_argument("--seed", type=int, default=2024)
    args = parser.parse_args()
    run_model(args.input, args.metadata, args.output_prefix, args.repeats, args.folds, args.seed)


if __name__ == "__main__":
    main()

