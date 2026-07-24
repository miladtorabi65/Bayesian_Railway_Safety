from pathlib import Path

import pandas as pd


# ============================================================
# 1. Paths
# ============================================================

DATA_PATH = "../data/raw/csidata_2006-2024.xlsx"

OUTPUT_DIR = Path("../data/processed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_PATH = (
    OUTPUT_DIR
    / "railway_country_year_2010_2024.csv"
)


# ============================================================
# 2. Project settings
# ============================================================

START_YEAR = 2010
END_YEAR = 2024
REFERENCE_YEAR = 2017

COUNTRIES = [
    "AT", "BE", "BG", "CH", "CZ", "DE", "DK",
    "EE", "EL", "ES", "FI", "FR", "HR", "HU",
    "IE", "IT", "LT", "LU", "LV", "NL", "NO",
    "PL", "PT", "RO", "SE", "SI", "SK",
]

COUNTRY_NAMES = {
    "AT": "Austria",
    "BE": "Belgium",
    "BG": "Bulgaria",
    "CH": "Switzerland",
    "CZ": "Czechia",
    "DE": "Germany",
    "DK": "Denmark",
    "EE": "Estonia",
    "EL": "Greece",
    "ES": "Spain",
    "FI": "Finland",
    "FR": "France",
    "HR": "Croatia",
    "HU": "Hungary",
    "IE": "Ireland",
    "IT": "Italy",
    "LT": "Lithuania",
    "LU": "Luxembourg",
    "LV": "Latvia",
    "NL": "Netherlands",
    "NO": "Norway",
    "PL": "Poland",
    "PT": "Portugal",
    "RO": "Romania",
    "SE": "Sweden",
    "SI": "Slovenia",
    "SK": "Slovakia",
}


# ============================================================
# 3. Load raw tables
# ============================================================

accident_raw = pd.read_excel(
    DATA_PATH,
    sheet_name="Table 0",
    header=1,
)

traffic_raw = pd.read_excel(
    DATA_PATH,
    sheet_name="Table 6",
    header=1,
)


# ============================================================
# 4. Fill indicator codes downward
# ============================================================

for df in [accident_raw, traffic_raw]:
    df["Code"] = df["Code"].ffill()
    df["CSI"] = df["CSI"].ffill()


# ============================================================
# 5. Extract relevant indicators
# ============================================================

accidents = accident_raw[
    (accident_raw["Code"] == "N00")
    & accident_raw["Year"].between(
        START_YEAR,
        END_YEAR,
    )
].copy()

train_km = traffic_raw[
    (traffic_raw["Code"] == "R01")
    & traffic_raw["Year"].between(
        START_YEAR,
        END_YEAR,
    )
].copy()


# ============================================================
# 6. Convert from wide to long format
# ============================================================

accidents_long = (
    accidents[
        ["Year"] + COUNTRIES
    ]
    .melt(
        id_vars="Year",
        var_name="country",
        value_name="accidents",
    )
    .rename(
        columns={"Year": "year"}
    )
)

train_km_long = (
    train_km[
        ["Year"] + COUNTRIES
    ]
    .melt(
        id_vars="Year",
        var_name="country",
        value_name="train_km",
    )
    .rename(
        columns={"Year": "year"}
    )
)


# ============================================================
# 7. Merge outcome and exposure
# ============================================================

data = accidents_long.merge(
    train_km_long,
    on=["country", "year"],
    how="inner",
    validate="one_to_one",
)

data["accidents"] = data["accidents"].astype(int)

# ============================================================
# 8. Add country names
# ============================================================

data["country_name"] = (
    data["country"]
    .map(COUNTRY_NAMES)
)


# ============================================================
# 9. Add analysis variables
# ============================================================

data["year_centered"] = (
    data["year"]
    - REFERENCE_YEAR
)

data[
    "accidents_per_million_train_km"
] = (
    data["accidents"]
    / data["train_km"]
)


# ============================================================
# 10. Add country IDs
# ============================================================

country_order = sorted(
    data["country"].unique()
)

country_to_id = {
    country: i + 1
    for i, country
    in enumerate(country_order)
}

data["country_id"] = (
    data["country"]
    .map(country_to_id)
)


# ============================================================
# 11. Sort data
# ============================================================

data = (
    data
    .sort_values(
        ["country_id", "year"]
    )
    .reset_index(drop=True)
)


# ============================================================
# 12. Validation
# ============================================================

assert len(data) == 405
assert data["country"].nunique() == 27
assert data["year"].nunique() == 15
assert data["year"].min() == 2010
assert data["year"].max() == 2024

assert not data.isna().any().any()

assert (
    data["train_km"] > 0
).all()

assert (
    data["accidents"] >= 0
).all()

assert (
    data["accidents"] % 1 == 0
).all()

assert (
    data.duplicated(
        ["country", "year"]
    ).sum()
    == 0
)

assert (
    data
    .groupby("country")
    .size()
    .eq(15)
    .all()
)


# ============================================================
# 13. Print summary
# ============================================================

print("=" * 50)
print("FINAL DATASET VALIDATION")
print("=" * 50)

print(
    f"Rows: "
    f"{len(data)}"
)

print(
    f"Countries: "
    f"{data['country'].nunique()}"
)

print(
    f"Years: "
    f"{data['year'].min()}-"
    f"{data['year'].max()}"
)

print(
    f"Missing values: "
    f"{data.isna().sum().sum()}"
)

print(
    f"Duplicate country-years: "
    f"{data.duplicated(['country', 'year']).sum()}"
)

print(
    f"Zero accident observations: "
    f"{(data['accidents'] == 0).sum()}"
)

print("\nFirst rows:")
print(data.head(10))


# ============================================================
# 14. Save frozen dataset
# ============================================================

data.to_csv(
    OUTPUT_PATH,
    index=False,
)

print(
    f"\nDataset saved to:\n"
    f"{OUTPUT_PATH}"
)