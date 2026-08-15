# =============================================================================
# PUBLICATION-GRADE TGARCH-X SENSITIVITY SUITE
# =============================================================================
# Purpose
#   Generate leakage-safe, one-trading-day-ahead sector variance forecasts for
#   the Nairobi Securities Exchange (or any similarly structured panel).
#
# Main design choices
#   1. Firm returns use lagged market-capitalization weights.
#   2. Every forecast uses only information available at the forecast origin.
#   3. The X regressor is lagged one trading observation.
#   4. Missing X values do not silently delete return dates.
#   5. TGARCH-X is compared with TGARCH without X, standard GARCH, GJR-GARCH,
#      EWMA, and rolling historical variance on common out-of-sample support.
#   6. Failed GARCH fits remain missing in the primary analysis. A clearly
#      labelled operational fallback is saved separately.
#   7. The common target and ratio-form QLIKE match the downstream Python code.
#   8. Pairwise loss tests use HAC inference, moving-block bootstrap intervals,
#      and Holm multiplicity adjustment.
#
# IMPORTANT
#   rugarch's fGARCH/TGARCH parameterization must be described exactly as
#   implemented by rugarch. Do not substitute a different threshold recursion
#   in the manuscript.
# =============================================================================

# =============================================================================
# 0. PACKAGES, CONFIGURATION, AND REPRODUCIBILITY
# =============================================================================
required_packages <- c(
  "dplyr",
  "tidyr",
  "lubridate",
  "rugarch",
  "progress",
  "zoo",
  "strucchange",
  "digest",
  "zip"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(rugarch)
  library(progress)
  library(zoo)
  library(strucchange)
  library(digest)
  library(zip)
})

SEED <- 42L
set.seed(SEED)

CODE_VERSION <- "tgarch_x_reviewer_aligned_sensitivity_v1_00_2016_2026"

PARENT_ENGINE_VERSION <- "tgarch_x_reviewer_aligned_v4_04_2016_2026"
PARENT_ENGINE_SHA256 <- "db0bf02c361c06a85aed817d299398ab0afd484e72c66584e0ba88198d2101c4"

# The authoritative PRIMARY R stage is frozen. This script performs only the
# reviewer-requested TGARCH-X sensitivity re-estimations. The primary
# 200-observation/Student-t/turnover-change-lag-1 specification is re-estimated
# once inside this suite solely as a within-run anchor for paired robustness
# comparisons; it does not replace the frozen primary outputs.
RUN_PROFILE <- "SENSITIVITY"
SENSITIVITY_FAMILIES <- c("WINDOW", "DISTRIBUTION", "X")

SAMPLE_START_DATE <- as.Date(Sys.getenv("NSE_SAMPLE_START", unset = "2016-08-01"))
SAMPLE_END_DATE <- as.Date(Sys.getenv("NSE_SAMPLE_END", unset = "2026-07-31"))
EXPECTED_N_FIRMS <- as.integer(Sys.getenv("NSE_EXPECTED_FIRMS", unset = "52"))
EXPECTED_N_SECTORS <- as.integer(Sys.getenv("NSE_EXPECTED_SECTORS", unset = "11"))

PRIMARY_WINDOW_SIZE <- 200L
WINDOW_SENSITIVITY <- c(125L, 200L, 250L)
PRIMARY_DISTRIBUTION <- "std"
DISTRIBUTION_SENSITIVITY <- c("norm", "std", "ged")
RETURN_SCALE <- 100
FORECAST_HORIZON <- 1L
EWMA_LAMBDA <- 0.94
TARGET_EPSILON <- 1e-8
FORECAST_EPSILON <- 1e-12
DIAGNOSTIC_LAG <- 10L
CAPPED_WEIGHT_MAX <- 0.25

EVALUATION_START_DATE <- as.Date(NA)
FINAL_TEST_START_DATE <- as.Date(NA)
SENSITIVITY_FINAL_TEST_START <- as.Date("2023-05-10")

MIN_VOLUME_COVERAGE <- 0.80
MIN_RANGE_COVERAGE <- 0.80
EXTREME_RETURN_THRESHOLD <- 0.50
ZERO_RETURN_TOLERANCE <- 1e-12
REQUIRE_CONTIGUOUS_SECTOR_CALENDAR <- FALSE

CALENDAR_GAP_JUSTIFICATION <- paste(
  "The source panel contains intermittent sector-date omissions and thin",
  "trading. No missing return is fabricated or zero-filled. The audit retains",
  "every detected gap, and affected forecasts must be described as forecasts",
  "for the next observed sector date rather than necessarily the next",
  "market-wide trading date."
)

PRIMARY_X_TRANSFORMATION <- "turnover_change"
PRIMARY_X_LAG <- 1L
X_SENSITIVITY_GRID <- data.frame(
  X_Transformation = c("turnover_change", "turnover_change", "sector_volume_change"),
  X_Lag = c(1L, 2L, 1L),
  stringsAsFactors = FALSE
)
ALLOW_X_IMPUTATION <- TRUE

DM_HAC_LAG <- NA_integer_
P_VALUE_ADJUSTMENT <- "holm"
RUN_STABILITY_ANALYSIS <- TRUE
STABILITY_ROLLING_WINDOW <- 60L
BOOTSTRAP_REPLICATIONS <- 5000L
BOOTSTRAP_BLOCK_LENGTH <- NA_integer_
MIN_TGARCH_X_VALID_COVERAGE <- 0.95

DATA_DIR <- Sys.getenv(
  "NSE_DATA_DIR",
  unset = "C:/..DATA From Office LapTop/PhD Research Paper/Data/Testing-Data Files/R-Data Import"
)
MASTER_FILENAME <- Sys.getenv(
  "NSE_MASTER_FILE",
  unset = "Final Master File_All variables 01082016-31072026.csv"
)

if (!dir.exists(DATA_DIR)) {
  stop("DATA_DIR does not exist: ", DATA_DIR,
       "\nSet NSE_DATA_DIR or edit DATA_DIR.")
}

PERIOD_TAG <- paste0(
  format(SAMPLE_START_DATE, "%Y%m%d"),
  "_",
  format(SAMPLE_END_DATE, "%Y%m%d")
)
OUTPUT_DIR <- file.path(
  DATA_DIR,
  paste0("TGARCH_X_Reviewer_Aligned_", PERIOD_TAG, "_SENSITIVITY")
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_csv_safe <- function(x, filename) {
  write.csv(x, file.path(OUTPUT_DIR, filename), row.names = FALSE, na = "")
}

cat("=============================================================================\n")
cat("REVIEWER-ALIGNED TGARCH-X AND ECONOMETRIC BENCHMARK FORECASTING\n")
cat("=============================================================================\n")
cat("Code version: ", CODE_VERSION, "\n", sep = "")
cat("Run profile: ", RUN_PROFILE, "\n", sep = "")
cat("Study period: ", SAMPLE_START_DATE, " to ", SAMPLE_END_DATE, "\n", sep = "")
cat("Data directory: ", normalizePath(DATA_DIR, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Output directory: ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Primary rolling window: ", PRIMARY_WINDOW_SIZE, " observations\n", sep = "")
cat("Primary innovation distribution: ", PRIMARY_DISTRIBUTION, "\n", sep = "")
cat("Primary X: ", PRIMARY_X_TRANSFORMATION, " at lag ", PRIMARY_X_LAG, "\n", sep = "")
cat("Allow X imputation: ", ALLOW_X_IMPUTATION, "\n", sep = "")
cat("=============================================================================\n\n")

# =============================================================================
# 1. LOAD AND STANDARDIZE THE MASTER FILE
# =============================================================================
csv_files <- list.files(
  DATA_DIR,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (nzchar(MASTER_FILENAME)) {
  master_file <- file.path(DATA_DIR, MASTER_FILENAME)
  if (!file.exists(master_file)) {
    stop("MASTER_FILENAME does not exist: ", master_file)
  }
} else {
  master_candidates <- csv_files[
    grepl("master", basename(csv_files), ignore.case = TRUE)
  ]
  
  if (length(master_candidates) == 0L) {
    stop(
      "No CSV filename containing 'Master' was found. ",
      "Set NSE_MASTER_FILE to the exact filename."
    )
  }
  
  if (length(master_candidates) > 1L) {
    stop(
      "Multiple master candidates were found: ",
      paste(basename(master_candidates), collapse = ", "),
      ". Set NSE_MASTER_FILE explicitly."
    )
  }
  
  master_file <- master_candidates[[1L]]
}

raw_df <- read.csv(
  master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Remove accidental UTF-8 BOMs and surrounding spaces from source headers.
# Repair blank, NA, and duplicate names before any dplyr operation. dplyr
# refuses to transform data frames that contain NA or empty column names.
original_header_names <- names(raw_df)
clean_header_names <- trimws(sub("^\ufeff", "", original_header_names))

blank_header <- is.na(clean_header_names) | clean_header_names == ""
if (any(blank_header)) {
  clean_header_names[blank_header] <- paste0(
    "Unnamed_Column_",
    which(blank_header)
  )
}

clean_header_names <- make.unique(clean_header_names, sep = "_Duplicate_")
names(raw_df) <- clean_header_names

header_name_audit <- data.frame(
  Column_Position = seq_along(clean_header_names),
  Original_Name = ifelse(
    is.na(original_header_names),
    "<NA>",
    original_header_names
  ),
  Repaired_Name = clean_header_names,
  Name_Repaired = ifelse(
    is.na(original_header_names),
    TRUE,
    original_header_names != clean_header_names
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(header_name_audit, "00_Header_Name_Repair_Audit.csv")

if (nrow(raw_df) == 0L) {
  stop("The master file contains no rows.")
}

canonical_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

resolve_column <- function(existing_names, aliases, label) {
  canonical_existing <- canonical_name(existing_names)
  alias_positions <- match(canonical_name(aliases), canonical_existing)
  alias_positions <- alias_positions[is.finite(alias_positions)]
  
  if (length(alias_positions) == 0L) {
    stop(
      "Could not identify the ", label, " column. Available columns: ",
      paste(existing_names, collapse = ", ")
    )
  }
  
  existing_names[[alias_positions[[1L]]]]
}

resolve_optional_column <- function(existing_names, aliases) {
  canonical_existing <- canonical_name(existing_names)
  alias_positions <- match(canonical_name(aliases), canonical_existing)
  alias_positions <- alias_positions[is.finite(alias_positions)]
  if (length(alias_positions) == 0L) return(NA_character_)
  existing_names[[alias_positions[[1L]]]]
}

date_source <- resolve_column(
  names(raw_df),
  c("Date", "Trading Date", "Trade Date"),
  "date"
)

sector_source <- resolve_column(
  names(raw_df),
  c("Sector", "Category", "Industry"),
  "sector"
)

stock_source <- resolve_column(
  names(raw_df),
  c("Stock", "Ticker", "Security", "Company", "Symbol"),
  "stock"
)

# Prefer an adjusted price when available.
close_source <- resolve_column(
  names(raw_df),
  c(
    "Adjusted Close", "Adjusted_Close", "Adj Close", "Adj_Close",
    "Close", "Close Price", "Closing Price"
  ),
  "close or adjusted-close"
)

volume_source <- resolve_column(
  names(raw_df),
  c("Trading Volume", "Trading_Volume", "Volume", "Share Volume"),
  "trading-volume"
)

market_cap_source <- resolve_column(
  names(raw_df),
  c(
    "Market Capitalization (KES)", "Market Capitalisation (KES)",
    "Market Capitalization", "Market_Capitalization",
    "Market Cap", "Market_Cap", "MarketCap"
  ),
  "market-capitalization"
)

high_source <- resolve_optional_column(
  names(raw_df),
  c("High", "High Price", "Daily High", "Highest Price")
)
low_source <- resolve_optional_column(
  names(raw_df),
  c("Low", "Low Price", "Daily Low", "Lowest Price")
)

raw_df$.__High__ <- if (!is.na(high_source)) raw_df[[high_source]] else NA_real_
raw_df$.__Low__ <- if (!is.na(low_source)) raw_df[[low_source]] else NA_real_

source_columns <- data.frame(
  Standard_Name = c(
    "Date", "Sector", "Stock", "Close", "Volume", "Market_Cap", "High", "Low"
  ),
  Source_Name = c(
    date_source, sector_source, stock_source, close_source, volume_source,
    market_cap_source,
    ifelse(is.na(high_source), "NOT_AVAILABLE", high_source),
    ifelse(is.na(low_source), "NOT_AVAILABLE", low_source)
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(source_columns, "00_Source_Column_Map.csv")

if (!grepl("adj", canonical_name(close_source), fixed = TRUE)) {
  warning(
    "The selected price column does not appear to be adjusted: ",
    close_source,
    ". Verify stock splits, rights issues, and dividends."
  )
}

parse_date_multi <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  
  if (inherits(x, c("POSIXct", "POSIXlt"))) {
    return(as.Date(x))
  }
  
  x_character <- trimws(as.character(x))
  parsed <- suppressWarnings(
    as.Date(
      parse_date_time(
        x_character,
        orders = c(
          "dmy", "ymd", "mdy",
          "d-b-Y", "d-B-Y",
          "Ymd", "dmY", "mdY"
        ),
        tz = "UTC",
        quiet = TRUE
      )
    )
  )
  
  # Handle plausible Excel serial dates that remain unresolved.
  unresolved <- is.na(parsed)
  numeric_candidate <- suppressWarnings(as.numeric(x_character))
  excel_like <- unresolved &
    is.finite(numeric_candidate) &
    numeric_candidate > 20000 &
    numeric_candidate < 80000
  
  if (any(excel_like)) {
    parsed[excel_like] <- as.Date(
      numeric_candidate[excel_like],
      origin = "1899-12-30"
    )
  }
  
  parsed
}

to_numeric_clean <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  
  x_character <- trimws(as.character(x))
  lower_character <- tolower(x_character)
  
  missing_token <- is.na(x) |
    x_character == "" |
    lower_character %in% c("na", "n/a", "null", "none", "-", "--")
  
  parenthesized_negative <- grepl(
    "^\\s*\\(.*\\)\\s*$",
    x_character,
    perl = TRUE
  )
  
  # Put the hyphen last in the bracket expression. This avoids TRE treating
  # it as an invalid character range, while still retaining negative signs.
  cleaned <- gsub(
    "[^0-9eE+.-]",
    "",
    x_character,
    perl = TRUE
  )
  
  invalid_cleaned <- cleaned %in% c(
    "", ".", "+", "-", "e", "E", "+.", "-."
  )
  cleaned[missing_token | invalid_cleaned] <- NA_character_
  
  values <- suppressWarnings(as.numeric(cleaned))
  
  negative_index <- parenthesized_negative & is.finite(values)
  values[negative_index] <- -abs(values[negative_index])
  
  values
}

df <- raw_df %>%
  transmute(
    Date = parse_date_multi(.data[[date_source]]),
    Sector = trimws(gsub("\\s+", " ", as.character(.data[[sector_source]]))),
    Stock = trimws(gsub("\\s+", " ", as.character(.data[[stock_source]]))),
    Close = to_numeric_clean(.data[[close_source]]),
    High = to_numeric_clean(.data[[".__High__"]]),
    Low = to_numeric_clean(.data[[".__Low__"]]),
    Trading_Volume = to_numeric_clean(.data[[volume_source]]),
    Market_Cap = to_numeric_clean(.data[[market_cap_source]])
  ) %>%
  mutate(
    Sector = case_when(
      tolower(Sector) == "commercial and services" ~ "Commercial and Services",
      TRUE ~ Sector
    )
  )

invalid_core_rows <- df %>%
  mutate(Source_Row = row_number()) %>%
  filter(
    is.na(Date) |
      is.na(Sector) |
      Sector == "" |
      is.na(Stock) |
      Stock == "" |
      !is.finite(Close) |
      Close <= 0 |
      !is.finite(Market_Cap) |
      Market_Cap <= 0
  )

write_csv_safe(invalid_core_rows, "01_Invalid_Core_Rows.csv")

df <- df %>%
  filter(
    !is.na(Date),
    Date >= SAMPLE_START_DATE,
    Date <= SAMPLE_END_DATE,
    !is.na(Sector),
    Sector != "",
    !is.na(Stock),
    Stock != "",
    is.finite(Close),
    Close > 0,
    is.finite(Market_Cap),
    Market_Cap > 0
  ) %>%
  arrange(Sector, Stock, Date)

duplicate_rows <- df %>%
  count(Sector, Stock, Date, name = "N") %>%
  filter(N > 1L)

write_csv_safe(duplicate_rows, "02_Duplicate_Keys.csv")

if (nrow(duplicate_rows) > 0L) {
  stop(
    "Duplicate Sector-Stock-Date observations remain after core cleaning. ",
    "Inspect 02_Duplicate_Keys.csv; do not average them silently."
  )
}

if (nrow(df) == 0L) {
  stop("No valid observations remain after core cleaning.")
}

# =============================================================================
# 2. FIRM RETURNS AND DATA AUDIT
# =============================================================================
max_run_length <- function(flag_vector) {
  flag_vector <- as.logical(flag_vector)
  if (length(flag_vector) == 0L || !any(flag_vector, na.rm = TRUE)) {
    return(0L)
  }
  runs <- rle(replace(flag_vector, is.na(flag_vector), FALSE))
  max(runs$lengths[runs$values])
}

trading_calendar <- data.frame(
  Date = sort(unique(df$Date)),
  Market_Trading_Index = seq_along(sort(unique(df$Date)))
)

df <- df %>%
  left_join(trading_calendar, by = "Date")

firm_panel_all <- df %>%
  group_by(Sector, Stock) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(
    Previous_Trading_Index = lag(Market_Trading_Index),
    Trading_Index_Gap = Market_Trading_Index - Previous_Trading_Index,
    Is_One_Trading_Day_Return = Trading_Index_Gap == 1L,
    Stock_Return = ifelse(
      Is_One_Trading_Day_Return,
      log(Close / lag(Close)),
      NA_real_
    ),
    Market_Cap_Lag1 = ifelse(
      Is_One_Trading_Day_Return,
      lag(Market_Cap),
      NA_real_
    ),
    Trading_Value = ifelse(
      is.finite(Trading_Volume) & Trading_Volume >= 0,
      Trading_Volume * Close,
      NA_real_
    ),
    Parkinson_Variance_Firm = ifelse(
      is.finite(High) & is.finite(Low) & High > 0 & Low > 0 & High >= Low,
      (log(High / Low)^2) / (4 * log(2)),
      NA_real_
    )
  ) %>%
  ungroup()

nonconsecutive_firm_observations <- firm_panel_all %>%
  filter(
    is.finite(Trading_Index_Gap),
    Trading_Index_Gap != 1L
  ) %>%
  select(
    Date,
    Sector,
    Stock,
    Trading_Index_Gap,
    Close,
    Market_Cap,
    Trading_Volume
  )

write_csv_safe(
  nonconsecutive_firm_observations,
  "03A_Nonconsecutive_Firm_Observations.csv"
)

firm_panel <- firm_panel_all %>%
  filter(
    is.finite(Stock_Return),
    is.finite(Market_Cap_Lag1),
    Market_Cap_Lag1 > 0
  )

if (nrow(firm_panel) == 0L) {
  stop("No valid one-trading-day firm returns could be constructed.")
}

firm_return_audit <- firm_panel %>%
  group_by(Sector, Stock) %>%
  summarise(
    First_Return_Date = min(Date),
    Last_Return_Date = max(Date),
    N_Returns = n(),
    Expected_Market_Return_Dates = max(nrow(trading_calendar) - 1L, 1L),
    Market_Date_Coverage = N_Returns / Expected_Market_Return_Dates,
    Meets_95pct_Coverage = Market_Date_Coverage >= 0.95,
    Zero_Return_Rate = mean(
      abs(Stock_Return) <= ZERO_RETURN_TOLERANCE,
      na.rm = TRUE
    ),
    Max_Zero_Return_Streak = max_run_length(
      abs(Stock_Return) <= ZERO_RETURN_TOLERANCE
    ),
    Mean_Return = mean(Stock_Return, na.rm = TRUE),
    SD_Return = sd(Stock_Return, na.rm = TRUE),
    Min_Return = min(Stock_Return, na.rm = TRUE),
    Max_Return = max(Stock_Return, na.rm = TRUE),
    N_Extreme_Returns = sum(
      abs(Stock_Return) > EXTREME_RETURN_THRESHOLD,
      na.rm = TRUE
    ),
    Missing_Volume_Rate = mean(
      !is.finite(Trading_Volume) | Trading_Volume < 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

extreme_firm_returns <- firm_panel %>%
  filter(abs(Stock_Return) > EXTREME_RETURN_THRESHOLD) %>%
  select(
    Date,
    Sector,
    Stock,
    Close,
    Stock_Return,
    Trading_Volume,
    Market_Cap,
    Market_Cap_Lag1
  ) %>%
  arrange(desc(abs(Stock_Return)))

write_csv_safe(firm_return_audit, "03_Firm_Return_Audit.csv")
write_csv_safe(extreme_firm_returns, "04_Extreme_Firm_Returns.csv")

# =============================================================================
# 3. SECTOR RETURNS, TURNOVER, AND THE ONE-DAY VARIANCE PROXY
# =============================================================================
sector_daily <- firm_panel %>%
  group_by(Date, Sector) %>%
  summarise(
    Weighted_Return_Sum = sum(Market_Cap_Lag1 * Stock_Return, na.rm = TRUE),
    Sector_Cap_Lag1 = sum(Market_Cap_Lag1, na.rm = TRUE),
    N_Stocks = n_distinct(Stock),
    N_Volume_Stocks = sum(is.finite(Trading_Volume) & Trading_Volume >= 0),
    Sector_Volume_Raw = if (
      sum(is.finite(Trading_Volume) & Trading_Volume >= 0) > 0L
    ) {
      sum(Trading_Volume[is.finite(Trading_Volume) & Trading_Volume >= 0],
          na.rm = TRUE)
    } else {
      NA_real_
    },
    Sector_Traded_Value_Raw = if (
      sum(is.finite(Trading_Value) & Trading_Value >= 0) > 0L
    ) {
      sum(Trading_Value[is.finite(Trading_Value) & Trading_Value >= 0],
          na.rm = TRUE)
    } else {
      NA_real_
    },
    N_Range_Stocks = sum(is.finite(Parkinson_Variance_Firm)),
    Range_Cap_Lag1 = sum(
      Market_Cap_Lag1[is.finite(Parkinson_Variance_Firm)],
      na.rm = TRUE
    ),
    Parkinson_Weighted_Sum = sum(
      Market_Cap_Lag1 * Parkinson_Variance_Firm,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Sector_Return = Weighted_Return_Sum / Sector_Cap_Lag1,
    Volume_Coverage = N_Volume_Stocks / N_Stocks,
    Range_Coverage = N_Range_Stocks / N_Stocks,
    Sector_Volume = ifelse(
      Volume_Coverage >= MIN_VOLUME_COVERAGE,
      Sector_Volume_Raw,
      NA_real_
    ),
    Sector_Traded_Value = ifelse(
      Volume_Coverage >= MIN_VOLUME_COVERAGE,
      Sector_Traded_Value_Raw,
      NA_real_
    ),
    Turnover_Ratio = Sector_Traded_Value / Sector_Cap_Lag1,
    Parkinson_Variance = ifelse(
      Range_Coverage >= MIN_RANGE_COVERAGE &
        is.finite(Range_Cap_Lag1) & Range_Cap_Lag1 > 0,
      Parkinson_Weighted_Sum / Range_Cap_Lag1,
      NA_real_
    )
  ) %>%
  filter(
    is.finite(Sector_Return),
    is.finite(Sector_Cap_Lag1),
    Sector_Cap_Lag1 > 0
  ) %>%
  arrange(Sector, Date)

if (anyDuplicated(sector_daily[, c("Sector", "Date")]) > 0L) {
  stop("Duplicate Sector-Date rows were created.")
}

sector_daily <- sector_daily %>%
  left_join(trading_calendar, by = "Date") %>%
  group_by(Sector) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(
    Sector_Trading_Index_Gap =
      Market_Trading_Index - lag(Market_Trading_Index)
  ) %>%
  ungroup()

sector_calendar_gaps <- sector_daily %>%
  filter(
    is.finite(Sector_Trading_Index_Gap),
    Sector_Trading_Index_Gap != 1L
  ) %>%
  select(
    Date,
    Sector,
    Sector_Trading_Index_Gap,
    N_Stocks,
    Volume_Coverage
  )

write_csv_safe(sector_calendar_gaps, "05A_Sector_Calendar_Gaps.csv")

sector_trading_sparsity_audit <- sector_daily %>%
  group_by(Sector) %>%
  summarise(
    N_Sector_Dates = n(),
    Min_Contributing_Stocks = min(N_Stocks, na.rm = TRUE),
    Median_Contributing_Stocks = median(N_Stocks, na.rm = TRUE),
    Mean_Contributing_Stocks = mean(N_Stocks, na.rm = TRUE),
    N_Dates_One_Stock = sum(N_Stocks == 1L, na.rm = TRUE),
    One_Stock_Date_Rate = mean(N_Stocks == 1L, na.rm = TRUE),
    N_Sector_Calendar_Gaps = sum(
      is.finite(Sector_Trading_Index_Gap) &
        Sector_Trading_Index_Gap != 1L,
      na.rm = TRUE
    ),
    Thin_Trading_Flag = (
      min(N_Stocks, na.rm = TRUE) <= 1L |
        mean(N_Stocks == 1L, na.rm = TRUE) >= 0.05
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(Thin_Trading_Flag), desc(One_Stock_Date_Rate), Sector)

write_csv_safe(
  sector_trading_sparsity_audit,
  "05B_Sector_Trading_Sparsity_Audit.csv"
)

sector_calendar_gap_summary <- if (nrow(sector_calendar_gaps) > 0L) {
  sector_calendar_gaps %>%
    group_by(Sector) %>%
    summarise(
      N_Gaps = n(),
      First_Gap_Date = min(Date),
      Last_Gap_Date = max(Date),
      Maximum_Trading_Index_Gap = max(Sector_Trading_Index_Gap),
      Mean_Trading_Index_Gap = mean(Sector_Trading_Index_Gap),
      .groups = "drop"
    ) %>%
    arrange(desc(N_Gaps), Sector)
} else {
  data.frame(
    Sector = character(0),
    N_Gaps = integer(0),
    First_Gap_Date = as.Date(character(0)),
    Last_Gap_Date = as.Date(character(0)),
    Maximum_Trading_Index_Gap = numeric(0),
    Mean_Trading_Index_Gap = numeric(0)
  )
}

write_csv_safe(
  sector_calendar_gap_summary,
  "05B_Sector_Calendar_Gap_Summary.csv"
)

if (nrow(sector_calendar_gaps) > 0L) {
  if (isTRUE(REQUIRE_CONTIGUOUS_SECTOR_CALENDAR)) {
    stop(
      "At least one sector series skips market trading dates. ",
      "Inspect 05A_Sector_Calendar_Gaps.csv and ",
      "05B_Sector_Calendar_Gap_Summary.csv. Repair the panel or explicitly ",
      "set REQUIRE_CONTIGUOUS_SECTOR_CALENDAR <- FALSE with justification."
    )
  }
  
  warning(
    "Continuing with ", nrow(sector_calendar_gaps),
    " nonconsecutive sector observations across ",
    n_distinct(sector_calendar_gaps$Sector), " sectors. ",
    CALENDAR_GAP_JUSTIFICATION,
    " Inspect outputs 05A and 05B."
  )
}

sector_calendar_audit <- sector_daily %>%
  group_by(Sector) %>%
  summarise(
    First_Date = min(Date),
    Last_Date = max(Date),
    N_Dates = n(),
    Max_Trading_Index_Gap = if (
      any(is.finite(Sector_Trading_Index_Gap))
    ) {
      max(
        Sector_Trading_Index_Gap[
          is.finite(Sector_Trading_Index_Gap)
        ]
      )
    } else {
      NA_real_
    },
    Mean_N_Stocks = mean(N_Stocks),
    Min_N_Stocks = min(N_Stocks),
    Max_N_Stocks = max(N_Stocks),
    Singleton_Date_Rate = mean(N_Stocks <= 1L),
    Mean_Volume_Coverage = mean(Volume_Coverage, na.rm = TRUE),
    Low_Volume_Coverage_Rate = mean(
      Volume_Coverage < MIN_VOLUME_COVERAGE,
      na.rm = TRUE
    ),
    Zero_Sector_Return_Rate = mean(
      abs(Sector_Return) <= ZERO_RETURN_TOLERANCE
    ),
    Max_Abs_Sector_Return = max(abs(Sector_Return)),
    .groups = "drop"
  )

write_csv_safe(sector_calendar_audit, "05_Sector_Calendar_Audit.csv")

singleton_sectors <- sector_calendar_audit %>%
  filter(Singleton_Date_Rate > 0)

if (nrow(singleton_sectors) > 0L) {
  warning(
    "At least one sector is represented by one stock on some dates. ",
    "Interpret those results as single-stock industry cases."
  )
}

common_actual <- sector_daily %>%
  transmute(
    Date,
    Sector,
    Sector_Return,
    Actual_Variance_Raw = Sector_Return^2,
    Actual_Variance = pmax(Sector_Return^2, TARGET_EPSILON),
    Parkinson_Variance = ifelse(
      is.finite(Parkinson_Variance) & Parkinson_Variance >= 0,
      pmax(Parkinson_Variance, TARGET_EPSILON),
      NA_real_
    ),
    N_Stocks,
    Sector_Cap_Lag1,
    Volume_Coverage,
    Range_Coverage,
    Code_Version = CODE_VERSION
  )

write_csv_safe(
  common_actual,
  "06_Common_Actual_Variance_1Day_Squared_Return.csv"
)

range_proxy_status <- data.frame(
  High_Column = ifelse(is.na(high_source), "NOT_AVAILABLE", high_source),
  Low_Column = ifelse(is.na(low_source), "NOT_AVAILABLE", low_source),
  N_Sector_Days_With_Parkinson = sum(
    is.finite(common_actual$Parkinson_Variance)
  ),
  Total_Sector_Days = nrow(common_actual),
  Coverage_Rate = mean(is.finite(common_actual$Parkinson_Variance)),
  Recommendation_Status = ifelse(
    mean(is.finite(common_actual$Parkinson_Variance)) >= MIN_RANGE_COVERAGE,
    "Available for range-based robustness",
    "Insufficient coverage; report data limitation"
  ),
  stringsAsFactors = FALSE
)
write_csv_safe(range_proxy_status, "06A_Alternative_Proxy_Status.csv")

# =============================================================================
# 4. LEAKAGE-SAFE X REGRESSOR
# =============================================================================
sector_model_data <- sector_daily %>%
  group_by(Sector) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(
    Log_Sector_Volume = ifelse(
      is.finite(Sector_Volume) & Sector_Volume >= 0,
      log1p(Sector_Volume),
      NA_real_
    ),
    Sector_Volume_Change = Log_Sector_Volume - lag(Log_Sector_Volume),
    Log_Turnover = ifelse(
      is.finite(Turnover_Ratio) & Turnover_Ratio >= 0,
      log1p(Turnover_Ratio),
      NA_real_
    ),
    Turnover_Change = Log_Turnover - lag(Log_Turnover),
    Turnover_Change_Lag1 = lag(Turnover_Change, 1L),
    Turnover_Change_Lag2 = lag(Turnover_Change, 2L),
    Sector_Volume_Change_Lag1 = lag(Sector_Volume_Change, 1L),
    Sector_Volume_Change_Lag2 = lag(Sector_Volume_Change, 2L),
    Actual_Variance_Raw = Sector_Return^2,
    Actual_Variance = pmax(Sector_Return^2, TARGET_EPSILON)
  ) %>%
  ungroup() %>%
  arrange(Sector, Date)

x_column_name <- function(transformation, lag_value) {
  if (identical(transformation, "turnover_change") && lag_value == 1L) {
    return("Turnover_Change_Lag1")
  }
  if (identical(transformation, "turnover_change") && lag_value == 2L) {
    return("Turnover_Change_Lag2")
  }
  if (identical(transformation, "sector_volume_change") && lag_value == 1L) {
    return("Sector_Volume_Change_Lag1")
  }
  if (identical(transformation, "sector_volume_change") && lag_value == 2L) {
    return("Sector_Volume_Change_Lag2")
  }
  stop("Unsupported X transformation/lag combination.")
}

x_audit_rows <- list()
for (x_index in seq_len(nrow(X_SENSITIVITY_GRID))) {
  transformation <- X_SENSITIVITY_GRID$X_Transformation[[x_index]]
  lag_value <- X_SENSITIVITY_GRID$X_Lag[[x_index]]
  column_name <- x_column_name(transformation, lag_value)
  audit_part <- sector_model_data %>%
    group_by(Sector) %>%
    summarise(
      X_Transformation = transformation,
      X_Lag = lag_value,
      X_Column = column_name,
      N_Observations = n(),
      N_X_Observed = sum(is.finite(.data[[column_name]])),
      X_Missing_Rate = mean(!is.finite(.data[[column_name]])),
      X_Mean = if (any(is.finite(.data[[column_name]]))) {
        mean(.data[[column_name]][is.finite(.data[[column_name]])])
      } else {
        NA_real_
      },
      X_SD = if (sum(is.finite(.data[[column_name]])) >= 2L) {
        sd(.data[[column_name]][is.finite(.data[[column_name]])])
      } else {
        NA_real_
      },
      X_Min = if (any(is.finite(.data[[column_name]]))) {
        min(.data[[column_name]][is.finite(.data[[column_name]])])
      } else {
        NA_real_
      },
      X_Max = if (any(is.finite(.data[[column_name]]))) {
        max(.data[[column_name]][is.finite(.data[[column_name]])])
      } else {
        NA_real_
      },
      .groups = "drop"
    )
  x_audit_rows[[x_index]] <- audit_part
}
x_audit <- bind_rows(x_audit_rows)
write_csv_safe(x_audit, "07_X_Regressor_Audit.csv")

eligible_sector_counts <- sector_model_data %>%
  count(Sector, name = "N_Observations") %>%
  mutate(
    Eligible_Primary = N_Observations > PRIMARY_WINDOW_SIZE,
    Eligible_Any_Sensitivity = N_Observations > min(WINDOW_SENSITIVITY)
  )

write_csv_safe(eligible_sector_counts, "08_Sector_Eligibility.csv")

eligible_sectors <- eligible_sector_counts %>%
  filter(Eligible_Any_Sensitivity) %>%
  pull(Sector)

if (length(eligible_sectors) == 0L) {
  stop("No sector has sufficient observations for any requested window.")
}

# Daily lagged-cap, equal, and exactly capped market weights.
cap_and_renormalize_weights <- function(raw_weights, cap_value) {
  w <- as.numeric(raw_weights)
  w[!is.finite(w) | w < 0] <- 0
  n <- length(w)
  if (n == 0L || sum(w) <= 0) return(rep(NA_real_, n))
  w <- w / sum(w)
  effective_cap <- max(cap_value, 1 / n)
  output <- rep(0, n)
  active <- rep(TRUE, n)
  remaining_mass <- 1
  
  while (any(active)) {
    active_weights <- w[active]
    if (sum(active_weights) <= 0) {
      output[active] <- remaining_mass / sum(active)
      break
    }
    proposal <- remaining_mass * active_weights / sum(active_weights)
    over <- proposal > effective_cap + 1e-12
    active_indices <- which(active)
    
    if (!any(over)) {
      output[active_indices] <- proposal
      break
    }
    
    capped_indices <- active_indices[over]
    output[capped_indices] <- effective_cap
    active[capped_indices] <- FALSE
    remaining_mass <- 1 - sum(output[!active])
    
    if (remaining_mass <= 1e-12) break
  }
  
  output / sum(output)
}

weight_panel <- sector_daily %>%
  group_by(Date) %>%
  arrange(Sector, .by_group = TRUE) %>%
  mutate(
    Market_Weight_Daily = Sector_Cap_Lag1 / sum(Sector_Cap_Lag1),
    Equal_Weight_Daily = 1 / n()
  ) %>%
  group_modify(
    ~{
      .x$Capped_Market_Weight_Daily <- cap_and_renormalize_weights(
        .x$Market_Weight_Daily,
        CAPPED_WEIGHT_MAX
      )
      .x
    }
  ) %>%
  ungroup() %>%
  select(
    Date, Sector, Sector_Cap_Lag1,
    Market_Weight_Daily, Equal_Weight_Daily,
    Capped_Market_Weight_Daily
  )

write_csv_safe(weight_panel, "08A_Daily_Weight_Schemes.csv")

weight_audit <- weight_panel %>%
  group_by(Date) %>%
  summarise(
    N_Sectors = n(),
    Market_Weight_Sum = sum(Market_Weight_Daily),
    Equal_Weight_Sum = sum(Equal_Weight_Daily),
    Capped_Weight_Sum = sum(Capped_Market_Weight_Daily),
    Maximum_Capped_Weight = max(Capped_Market_Weight_Daily),
    .groups = "drop"
  )
write_csv_safe(weight_audit, "08B_Daily_Weight_Audit.csv")

if (
  any(abs(weight_audit$Market_Weight_Sum - 1) > 1e-10) ||
  any(abs(weight_audit$Equal_Weight_Sum - 1) > 1e-10) ||
  any(abs(weight_audit$Capped_Weight_Sum - 1) > 1e-10)
) {
  stop("At least one daily weight scheme does not sum to one.")
}

# =============================================================================
# 5. FORECASTING HELPERS AND SPECIFICATION GRID
# =============================================================================
safe_variance <- function(x, epsilon = FORECAST_EPSILON) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  value <- suppressWarnings(var(x))
  if (!is.finite(value) || value <= epsilon) {
    value <- mean(x^2)
  }
  if (!is.finite(value) || value <= epsilon) value <- epsilon
  max(value, epsilon)
}

ewma_forecast <- function(x, lambda = EWMA_LAMBDA) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  variance_state <- safe_variance(x)
  for (i in seq_along(x)) {
    variance_state <- lambda * variance_state + (1 - lambda) * x[[i]]^2
  }
  max(variance_state, FORECAST_EPSILON)
}

standardize_training_vector <- function(x) {
  x <- as.numeric(x)
  center <- mean(x)
  scale_value <- sd(x)
  if (!is.finite(center)) stop("Training regressor mean is non-finite.")
  if (!is.finite(scale_value) || scale_value < 1e-12) scale_value <- 1
  list(
    values = (x - center) / scale_value,
    center = center,
    scale = scale_value
  )
}

prepare_x_window <- function(training_x, next_x) {
  training_x <- as.numeric(training_x)
  next_x <- as.numeric(next_x)
  missing_training <- !is.finite(training_x)
  missing_next <- length(next_x) != 1L || !is.finite(next_x)
  
  if (!ALLOW_X_IMPUTATION && (any(missing_training) || missing_next)) {
    return(list(
      valid = FALSE,
      reason = "missing_x_in_training_or_forecast_origin"
    ))
  }
  
  if (ALLOW_X_IMPUTATION) {
    observed_training <- training_x[is.finite(training_x)]
    if (length(observed_training) == 0L) {
      return(list(valid = FALSE, reason = "all_training_x_missing"))
    }
    imputation_value <- median(observed_training)
    training_imputed <- training_x
    training_imputed[missing_training] <- imputation_value
    next_imputed <- if (missing_next) imputation_value else next_x
    scaling <- standardize_training_vector(training_imputed)
    
    x_train <- cbind(
      X_z = scaling$values,
      X_Missing = as.numeric(missing_training)
    )
    x_next <- matrix(
      c(
        (next_imputed - scaling$center) / scaling$scale,
        as.numeric(missing_next)
      ),
      nrow = 1L
    )
    colnames(x_next) <- colnames(x_train)
    
    return(list(
      valid = TRUE,
      reason = NA_character_,
      x_train = x_train,
      x_next = x_next,
      center = scaling$center,
      scale = scaling$scale,
      imputation_value = imputation_value,
      used_imputation = any(missing_training) || missing_next,
      n_missing_training = sum(missing_training),
      missing_forecast_origin = missing_next,
      training_imputation_rate = mean(missing_training)
    ))
  }
  
  scaling <- standardize_training_vector(training_x)
  x_train <- matrix(
    scaling$values,
    ncol = 1L,
    dimnames = list(NULL, "X_z")
  )
  x_next <- matrix(
    (next_x - scaling$center) / scaling$scale,
    nrow = 1L,
    dimnames = list(NULL, "X_z")
  )
  
  list(
    valid = TRUE,
    reason = NA_character_,
    x_train = x_train,
    x_next = x_next,
    center = scaling$center,
    scale = scaling$scale,
    imputation_value = NA_real_,
    used_imputation = FALSE,
    n_missing_training = 0L,
    missing_forecast_origin = FALSE,
    training_imputation_rate = 0
  )
}

sanitize_id <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

make_specification_id <- function(
    model_family,
    window_size,
    distribution,
    x_transformation,
    x_lag
) {
  parts <- c(
    model_family,
    paste0("w", window_size),
    ifelse(is.na(distribution), "no_dist", distribution),
    ifelse(
      is.na(x_transformation),
      "no_x",
      paste0(x_transformation, "_lag", x_lag)
    )
  )
  sanitize_id(paste(parts, collapse = "__"))
}

model_display_label <- function(model_family) {
  switch(
    model_family,
    "RollingVariance" = "Historical variance",
    "EWMA" = "EWMA",
    "sGARCH" = "GARCH(1,1)",
    "gjrGARCH" = "GJR-GARCH(1,1)",
    "fGARCH-TGARCH" = "TGARCH(1,1)",
    "fGARCH-TGARCH-X" = "TGARCH-X(1,1)",
    model_family
  )
}

build_specification_grid <- function() {
  primary <- data.frame(
    Model_Family = c(
      "RollingVariance",
      "EWMA",
      "sGARCH",
      "gjrGARCH",
      "fGARCH-TGARCH",
      "fGARCH-TGARCH-X"
    ),
    Window_Size = rep(PRIMARY_WINDOW_SIZE, 6L),
    Distribution = c(
      NA_character_,
      NA_character_,
      rep(PRIMARY_DISTRIBUTION, 4L)
    ),
    X_Transformation = c(
      rep(NA_character_, 5L),
      PRIMARY_X_TRANSFORMATION
    ),
    X_Lag = c(rep(NA_integer_, 5L), PRIMARY_X_LAG),
    Analysis_Set = "Primary",
    stringsAsFactors = FALSE
  )
  
  window_grid <- data.frame(
    Model_Family = "fGARCH-TGARCH-X",
    Window_Size = WINDOW_SENSITIVITY,
    Distribution = PRIMARY_DISTRIBUTION,
    X_Transformation = PRIMARY_X_TRANSFORMATION,
    X_Lag = PRIMARY_X_LAG,
    Analysis_Set = "Window_Sensitivity",
    stringsAsFactors = FALSE
  )
  
  distribution_grid <- data.frame(
    Model_Family = "fGARCH-TGARCH-X",
    Window_Size = PRIMARY_WINDOW_SIZE,
    Distribution = DISTRIBUTION_SENSITIVITY,
    X_Transformation = PRIMARY_X_TRANSFORMATION,
    X_Lag = PRIMARY_X_LAG,
    Analysis_Set = "Distribution_Sensitivity",
    stringsAsFactors = FALSE
  )
  
  x_grid <- X_SENSITIVITY_GRID %>%
    mutate(
      Model_Family = "fGARCH-TGARCH-X",
      Window_Size = PRIMARY_WINDOW_SIZE,
      Distribution = PRIMARY_DISTRIBUTION,
      Analysis_Set = "X_Sensitivity"
    ) %>%
    select(
      Model_Family,
      Window_Size,
      Distribution,
      X_Transformation,
      X_Lag,
      Analysis_Set
    )
  
  # Run all three sensitivity families together. The primary-anchor
  # TGARCH-X specification appears in each family but is deduplicated below
  # by Specification_ID so it is estimated only once.
  selected <- bind_rows(
    window_grid,
    distribution_grid,
    x_grid
  )
  
  selected %>%
    mutate(
      Specification_ID = mapply(
        make_specification_id,
        Model_Family,
        Window_Size,
        Distribution,
        X_Transformation,
        X_Lag,
        USE.NAMES = FALSE
      )
    ) %>%
    group_by(
      Specification_ID,
      Model_Family,
      Window_Size,
      Distribution,
      X_Transformation,
      X_Lag
    ) %>%
    summarise(
      Analysis_Set = paste(sort(unique(Analysis_Set)), collapse = "|"),
      .groups = "drop"
    ) %>%
    mutate(
      Model_Label = vapply(
        Model_Family,
        model_display_label,
        character(1)
      ),
      Uses_X = Model_Family == "fGARCH-TGARCH-X",
      Reestimation_Type = ifelse(
        grepl("Primary", Analysis_Set, fixed = TRUE),
        "Primary full estimation",
        "Sensitivity full re-estimation"
      ),
      Primary_Flag =
        Model_Family %in% primary$Model_Family &
        Window_Size == PRIMARY_WINDOW_SIZE &
        (
          is.na(Distribution) |
            Distribution == PRIMARY_DISTRIBUTION
        ) &
        (
          !Uses_X |
            (
              X_Transformation == PRIMARY_X_TRANSFORMATION &
                X_Lag == PRIMARY_X_LAG
            )
        )
    ) %>%
    arrange(Model_Family, Window_Size, Distribution, X_Transformation, X_Lag)
}

specification_grid <- build_specification_grid()
write_csv_safe(specification_grid, "09_Specification_Grid.csv")

make_garch_spec <- function(
    model_family,
    distribution,
    x_train = NULL
) {
  if (!distribution %in% c("norm", "std", "ged")) {
    stop("Unsupported innovation distribution: ", distribution)
  }
  
  variance_model <- switch(
    model_family,
    "sGARCH" = list(
      model = "sGARCH",
      garchOrder = c(1L, 1L),
      variance.targeting = FALSE
    ),
    "gjrGARCH" = list(
      model = "gjrGARCH",
      garchOrder = c(1L, 1L),
      variance.targeting = FALSE
    ),
    "fGARCH-TGARCH" = list(
      model = "fGARCH",
      submodel = "TGARCH",
      garchOrder = c(1L, 1L),
      variance.targeting = FALSE
    ),
    "fGARCH-TGARCH-X" = list(
      model = "fGARCH",
      submodel = "TGARCH",
      garchOrder = c(1L, 1L),
      external.regressors = x_train,
      variance.targeting = FALSE
    ),
    stop("Unknown GARCH model family: ", model_family)
  )
  
  ugarchspec(
    mean.model = list(
      armaOrder = c(0L, 0L),
      include.mean = TRUE
    ),
    variance.model = variance_model,
    distribution.model = distribution
  )
}

capture_ugarch_fit <- function(specification, data_vector) {
  warning_messages <- character(0)
  error_message <- NA_character_
  
  fit_object <- tryCatch(
    withCallingHandlers(
      ugarchfit(
        spec = specification,
        data = data_vector,
        solver = "hybrid",
        solver.control = list(trace = 0),
        fit.control = list(stationarity = 1, scale = 0)
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  
  list(
    fit = fit_object,
    warning = if (length(warning_messages) > 0L) {
      paste(unique(warning_messages), collapse = " | ")
    } else {
      NA_character_
    },
    error = error_message
  )
}

empty_fit_result <- function(
    failure_reason,
    convergence_code = NA_integer_,
    warning = NA_character_,
    error = NA_character_,
    fit = NULL
) {
  list(
    valid = FALSE,
    forecast = NA_real_,
    fit = fit,
    convergence_code = convergence_code,
    warning = warning,
    error = error,
    failure_reason = failure_reason,
    persistence = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    Log_Likelihood = NA_real_
  )
}

fit_and_forecast_garch <- function(
    model_family,
    distribution,
    y_training,
    x_prepared = NULL
) {
  x_train <- NULL
  x_next <- NULL
  
  if (identical(model_family, "fGARCH-TGARCH-X")) {
    if (is.null(x_prepared) || !isTRUE(x_prepared$valid)) {
      return(
        empty_fit_result(
          if (is.null(x_prepared)) {
            "x_preparation_not_supplied"
          } else {
            x_prepared$reason
          }
        )
      )
    }
    x_train <- x_prepared$x_train
    x_next <- x_prepared$x_next
  }
  
  specification <- tryCatch(
    make_garch_spec(model_family, distribution, x_train),
    error = function(e) e
  )
  if (inherits(specification, "error")) {
    return(empty_fit_result("specification_error", error = conditionMessage(specification)))
  }
  
  fit_capture <- capture_ugarch_fit(
    specification,
    y_training * RETURN_SCALE
  )
  fit_object <- fit_capture$fit
  
  if (is.null(fit_object)) {
    return(
      empty_fit_result(
        "fit_error",
        warning = fit_capture$warning,
        error = fit_capture$error
      )
    )
  }
  
  convergence_code <- tryCatch(
    as.integer(convergence(fit_object)),
    error = function(e) NA_integer_
  )
  
  if (!is.finite(convergence_code) || convergence_code != 0L) {
    return(
      empty_fit_result(
        "fit_nonconvergence",
        convergence_code = convergence_code,
        warning = fit_capture$warning,
        error = fit_capture$error,
        fit = fit_object
      )
    )
  }
  
  forecast_arguments <- list(fitORspec = fit_object, n.ahead = 1L)
  if (identical(model_family, "fGARCH-TGARCH-X")) {
    forecast_arguments$external.forecasts <- list(vregfor = x_next)
  }
  
  forecast_object <- tryCatch(
    do.call(ugarchforecast, forecast_arguments),
    error = function(e) e
  )
  if (inherits(forecast_object, "error")) {
    return(
      empty_fit_result(
        "forecast_error",
        convergence_code = convergence_code,
        warning = fit_capture$warning,
        error = conditionMessage(forecast_object),
        fit = fit_object
      )
    )
  }
  
  sigma_scaled <- tryCatch(
    as.numeric(sigma(forecast_object)[1L]),
    error = function(e) NA_real_
  )
  forecast_variance <- (sigma_scaled / RETURN_SCALE)^2
  
  if (
    length(forecast_variance) != 1L ||
    !is.finite(forecast_variance) ||
    forecast_variance <= 0
  ) {
    return(
      empty_fit_result(
        "nonfinite_or_nonpositive_forecast",
        convergence_code = convergence_code,
        warning = fit_capture$warning,
        error = fit_capture$error,
        fit = fit_object
      )
    )
  }
  
  information_criteria <- tryCatch(
    as.numeric(rugarch::infocriteria(fit_object)),
    error = function(e) rep(NA_real_, 4L)
  )
  
  list(
    valid = TRUE,
    forecast = max(forecast_variance, FORECAST_EPSILON),
    fit = fit_object,
    convergence_code = convergence_code,
    warning = fit_capture$warning,
    error = fit_capture$error,
    failure_reason = NA_character_,
    persistence = tryCatch(
      as.numeric(rugarch::persistence(fit_object)),
      error = function(e) NA_real_
    ),
    AIC = if (length(information_criteria) >= 1L) {
      information_criteria[[1L]]
    } else {
      NA_real_
    },
    BIC = if (length(information_criteria) >= 2L) {
      information_criteria[[2L]]
    } else {
      NA_real_
    },
    Log_Likelihood = tryCatch(
      as.numeric(likelihood(fit_object)),
      error = function(e) NA_real_
    )
  )
}

evaluation_phase <- function(date_value) {
  if (!is.na(FINAL_TEST_START_DATE) && date_value >= FINAL_TEST_START_DATE) {
    return("final_test")
  }
  if (!is.na(EVALUATION_START_DATE) && date_value >= EVALUATION_START_DATE) {
    return("development")
  }
  if (is.na(EVALUATION_START_DATE) && is.na(FINAL_TEST_START_DATE)) {
    return("evaluation")
  }
  "pre_evaluation"
}

extract_parameter_rows <- function(
    fit_object,
    sector_name,
    forecast_date,
    specification_row
) {
  if (is.null(fit_object)) return(NULL)
  estimates <- tryCatch(coef(fit_object), error = function(e) NULL)
  if (is.null(estimates) || length(estimates) == 0L) return(NULL)
  
  data.frame(
    Sector = sector_name,
    Forecast_Date = forecast_date,
    Specification_ID = specification_row$Specification_ID,
    Model_Family = specification_row$Model_Family,
    Window_Size = specification_row$Window_Size,
    Distribution = specification_row$Distribution,
    X_Transformation = specification_row$X_Transformation,
    X_Lag = specification_row$X_Lag,
    Parameter = names(estimates),
    Estimate = as.numeric(estimates),
    Code_Version = CODE_VERSION,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 6. ROLLING FORECAST ENGINE
# =============================================================================
sector_row_counts <- sector_model_data %>%
  count(Sector, name = "N")

total_origins <- 0L
for (spec_index in seq_len(nrow(specification_grid))) {
  window_size <- specification_grid$Window_Size[[spec_index]]
  total_origins <- total_origins + sum(
    pmax(sector_row_counts$N - window_size, 0L)
  )
}

if (total_origins <= 0L) {
  stop("No feasible rolling forecast origins were found.")
}

progress <- progress_bar$new(
  total = total_origins,
  clear = FALSE,
  width = 100,
  format = paste0(
    "[:bar] :percent | origins :current/:total | ",
    "valid GARCH :valid | failed :failed | ETA :eta"
  )
)

results <- vector("list", total_origins)
parameter_rows <- list()
last_successful_fit <- list()
last_successful_metadata <- list()

result_index <- 0L
parameter_index <- 0L
valid_garch_count <- 0L
failed_garch_count <- 0L
start_time <- Sys.time()

for (spec_index in seq_len(nrow(specification_grid))) {
  specification_row <- specification_grid[spec_index, , drop = FALSE]
  model_family <- specification_row$Model_Family[[1L]]
  window_size <- specification_row$Window_Size[[1L]]
  distribution <- specification_row$Distribution[[1L]]
  uses_x <- isTRUE(specification_row$Uses_X[[1L]])
  
  x_column <- if (uses_x) {
    x_column_name(
      specification_row$X_Transformation[[1L]],
      specification_row$X_Lag[[1L]]
    )
  } else {
    NA_character_
  }
  
  for (sector_name in eligible_sectors) {
    sector_data <- sector_model_data %>%
      filter(Sector == sector_name) %>%
      arrange(Date)
    
    n_sector <- nrow(sector_data)
    n_origins <- n_sector - window_size
    if (n_origins <= 0L) next
    
    for (origin_index in seq_len(n_origins)) {
      training_indices <- origin_index:(origin_index + window_size - 1L)
      forecast_index <- origin_index + window_size
      training_data <- sector_data[training_indices, , drop = FALSE]
      forecast_row <- sector_data[forecast_index, , drop = FALSE]
      
      if (
        !is.na(EVALUATION_START_DATE) &&
        forecast_row$Date < EVALUATION_START_DATE
      ) {
        progress$tick(
          tokens = list(valid = valid_garch_count, failed = failed_garch_count)
        )
        next
      }
      
      y_training <- training_data$Sector_Return
      historical_variance <- safe_variance(y_training)
      x_prepared <- if (uses_x) {
        prepare_x_window(
          training_data[[x_column]],
          forecast_row[[x_column]]
        )
      } else {
        NULL
      }
      
      is_simple <- model_family %in% c("RollingVariance", "EWMA")
      fitted <- NULL
      
      if (identical(model_family, "RollingVariance")) {
        forecast_value <- historical_variance
        forecast_valid <- is.finite(forecast_value) && forecast_value > 0
        failure_reason <- NA_character_
      } else if (identical(model_family, "EWMA")) {
        forecast_value <- ewma_forecast(y_training)
        forecast_valid <- is.finite(forecast_value) && forecast_value > 0
        failure_reason <- if (forecast_valid) NA_character_ else "ewma_failure"
      } else {
        fitted <- fit_and_forecast_garch(
          model_family = model_family,
          distribution = distribution,
          y_training = y_training,
          x_prepared = x_prepared
        )
        forecast_value <- if (isTRUE(fitted$valid)) fitted$forecast else NA_real_
        forecast_valid <- isTRUE(fitted$valid)
        failure_reason <- fitted$failure_reason
        
        if (forecast_valid) {
          valid_garch_count <- valid_garch_count + 1L
          fit_key <- paste(
            sector_name,
            specification_row$Specification_ID[[1L]],
            sep = "|||"
          )
          last_successful_fit[[fit_key]] <- fitted$fit
          last_successful_metadata[[fit_key]] <- data.frame(
            Sector = sector_name,
            Specification_ID = specification_row$Specification_ID[[1L]],
            Model_Family = model_family,
            Window_Start_Date = min(training_data$Date),
            Window_End_Date = max(training_data$Date),
            Forecast_Date = forecast_row$Date,
            stringsAsFactors = FALSE
          )
          
          if (RUN_PROFILE %in% c("PRIMARY", "SENSITIVITY")) {
            parameter_output <- extract_parameter_rows(
              fitted$fit,
              sector_name,
              forecast_row$Date,
              specification_row
            )
            if (!is.null(parameter_output)) {
              parameter_index <- parameter_index + 1L
              parameter_rows[[parameter_index]] <- parameter_output
            }
          }
        } else {
          failed_garch_count <- failed_garch_count + 1L
        }
      }
      
      operational_forecast <- if (forecast_valid) {
        forecast_value
      } else {
        historical_variance
      }
      
      result_index <- result_index + 1L
      results[[result_index]] <- data.frame(
        Date = forecast_row$Date,
        Sector = sector_name,
        Evaluation_Phase = evaluation_phase(forecast_row$Date),
        Window_Start_Date = min(training_data$Date),
        Window_End_Date = max(training_data$Date),
        Window_Observations = nrow(training_data),
        Specification_ID = specification_row$Specification_ID[[1L]],
        Analysis_Set = specification_row$Analysis_Set[[1L]],
        Reestimation_Type = specification_row$Reestimation_Type[[1L]],
        Model_Family = model_family,
        Model_Label = specification_row$Model_Label[[1L]],
        Distribution = distribution,
        X_Transformation = specification_row$X_Transformation[[1L]],
        X_Lag = specification_row$X_Lag[[1L]],
        Primary_Flag = specification_row$Primary_Flag[[1L]],
        Forecast_Variance = forecast_value,
        Operational_Forecast_Variance = operational_forecast,
        Fallback_Variance = historical_variance,
        Actual_Variance_Raw = forecast_row$Actual_Variance_Raw,
        Actual_Variance = forecast_row$Actual_Variance,
        Parkinson_Variance = forecast_row$Parkinson_Variance,
        Sector_Return = forecast_row$Sector_Return,
        Sector_Cap_Lag1 = forecast_row$Sector_Cap_Lag1,
        N_Stocks = forecast_row$N_Stocks,
        Volume_Coverage = forecast_row$Volume_Coverage,
        Range_Coverage = forecast_row$Range_Coverage,
        X_Value = if (uses_x) forecast_row[[x_column]] else NA_real_,
        X_Training_Mean = if (uses_x && isTRUE(x_prepared$valid)) {
          x_prepared$center
        } else {
          NA_real_
        },
        X_Training_SD = if (uses_x && isTRUE(x_prepared$valid)) {
          x_prepared$scale
        } else {
          NA_real_
        },
        X_Imputation_Value = if (uses_x && isTRUE(x_prepared$valid)) {
          x_prepared$imputation_value
        } else {
          NA_real_
        },
        X_Imputation_Used = if (uses_x && isTRUE(x_prepared$valid)) {
          x_prepared$used_imputation
        } else {
          FALSE
        },
        X_Training_Missing_Count = if (uses_x && isTRUE(x_prepared$valid)) {
          x_prepared$n_missing_training
        } else {
          NA_integer_
        },
        X_Training_Imputation_Rate = if (uses_x && isTRUE(x_prepared$valid)) {
          x_prepared$training_imputation_rate
        } else {
          NA_real_
        },
        X_Forecast_Origin_Imputed = if (uses_x && isTRUE(x_prepared$valid)) {
          x_prepared$missing_forecast_origin
        } else {
          FALSE
        },
        Convergence_Code = if (is_simple) NA_integer_ else fitted$convergence_code,
        Converged = if (is_simple) NA else forecast_valid,
        Forecast_Valid = forecast_valid,
        Operational_Fallback_Used = !forecast_valid,
        Failure_Reason = failure_reason,
        Fit_Warning = if (is_simple) NA_character_ else fitted$warning,
        Fit_Error = if (is_simple) NA_character_ else fitted$error,
        Persistence = if (is_simple) NA_real_ else fitted$persistence,
        AIC = if (is_simple) NA_real_ else fitted$AIC,
        BIC = if (is_simple) NA_real_ else fitted$BIC,
        Log_Likelihood = if (is_simple) NA_real_ else fitted$Log_Likelihood,
        Forecast_Horizon = FORECAST_HORIZON,
        Code_Version = CODE_VERSION,
        stringsAsFactors = FALSE
      )
      
      progress$tick(
        tokens = list(valid = valid_garch_count, failed = failed_garch_count)
      )
    }
  }
}

end_time <- Sys.time()

if (result_index == 0L) {
  stop("The rolling engine produced no forecasts.")
}

results_df <- bind_rows(results[seq_len(result_index)]) %>%
  arrange(Specification_ID, Sector, Date)

duplicate_forecasts <- results_df %>%
  count(Specification_ID, Sector, Date, name = "N") %>%
  filter(N > 1L)

if (nrow(duplicate_forecasts) > 0L) {
  write_csv_safe(duplicate_forecasts, "10A_Duplicate_Forecasts.csv")
  stop("Duplicate Specification-Sector-Date forecasts were generated.")
}

parameter_df <- if (length(parameter_rows) > 0L) {
  bind_rows(parameter_rows)
} else {
  data.frame(Message = "Rolling parameters were not saved for this profile.")
}

write_csv_safe(results_df, "10_All_Rolling_Forecasts_Long.csv")
write_csv_safe(parameter_df, "11_Rolling_GARCH_Parameters_Long.csv")

# =============================================================================
# 7. FIT COVERAGE, FAILURES, AND FIT QUALITY
# =============================================================================
fit_summary <- results_df %>%
  group_by(
    Specification_ID,
    Model_Family,
    Window_Observations,
    Distribution,
    X_Transformation,
    X_Lag,
    Sector,
    Evaluation_Phase
  ) %>%
  summarise(
    N_Origins = n(),
    N_Valid = sum(Forecast_Valid, na.rm = TRUE),
    Valid_Rate = mean(Forecast_Valid, na.rm = TRUE),
    N_Operational_Fallback = sum(Operational_Fallback_Used, na.rm = TRUE),
    Operational_Fallback_Rate = mean(
      Operational_Fallback_Used,
      na.rm = TRUE
    ),
    Mean_Persistence = mean(Persistence, na.rm = TRUE),
    Maximum_Persistence = if (any(is.finite(Persistence))) {
      max(Persistence, na.rm = TRUE)
    } else {
      NA_real_
    },
    Mean_AIC = mean(AIC, na.rm = TRUE),
    Mean_BIC = mean(BIC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_Persistence = ifelse(
      is.nan(Mean_Persistence),
      NA_real_,
      Mean_Persistence
    ),
    Mean_AIC = ifelse(is.nan(Mean_AIC), NA_real_, Mean_AIC),
    Mean_BIC = ifelse(is.nan(Mean_BIC), NA_real_, Mean_BIC)
  )

failure_summary <- results_df %>%
  filter(!Forecast_Valid) %>%
  mutate(
    Failure_Reason = ifelse(
      is.na(Failure_Reason),
      "unspecified",
      Failure_Reason
    )
  ) %>%
  count(
    Specification_ID,
    Model_Family,
    Sector,
    Evaluation_Phase,
    Failure_Reason,
    name = "N"
  ) %>%
  arrange(Specification_ID, Sector, desc(N))

write_csv_safe(fit_summary, "12_Fit_Coverage_and_Quality.csv")
write_csv_safe(failure_summary, "13_Failure_Reasons.csv")

# Reviewer-requested transparency for missing external-regressor values.
# The median used for imputation is calculated from the current rolling
# training window only. No future observation enters the imputation value.
x_imputation_audit <- results_df %>%
  filter(Model_Family == "fGARCH-TGARCH-X") %>%
  group_by(
    Specification_ID,
    Sector,
    Evaluation_Phase
  ) %>%
  summarise(
    N_Origins = n(),
    N_Valid_Forecasts = sum(Forecast_Valid, na.rm = TRUE),
    Valid_Coverage = mean(Forecast_Valid, na.rm = TRUE),
    N_Windows_With_Any_Imputation = sum(X_Imputation_Used, na.rm = TRUE),
    Window_Imputation_Rate = mean(X_Imputation_Used, na.rm = TRUE),
    Total_Training_X_Values_Imputed = sum(
      X_Training_Missing_Count,
      na.rm = TRUE
    ),
    Mean_Training_Imputation_Rate = mean(
      X_Training_Imputation_Rate,
      na.rm = TRUE
    ),
    N_Forecast_Origins_Imputed = sum(
      X_Forecast_Origin_Imputed,
      na.rm = TRUE
    ),
    Forecast_Origin_Imputation_Rate = mean(
      X_Forecast_Origin_Imputed,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_Training_Imputation_Rate = ifelse(
      is.nan(Mean_Training_Imputation_Rate),
      NA_real_,
      Mean_Training_Imputation_Rate
    )
  ) %>%
  arrange(Specification_ID, Sector, Evaluation_Phase)

write_csv_safe(x_imputation_audit, "13A_TGARCH_X_Imputation_Audit.csv")

primary_tgarch_x_id <- specification_grid %>%
  filter(
    Model_Family == "fGARCH-TGARCH-X",
    Window_Size == PRIMARY_WINDOW_SIZE,
    Distribution == PRIMARY_DISTRIBUTION,
    X_Transformation == PRIMARY_X_TRANSFORMATION,
    X_Lag == PRIMARY_X_LAG
  ) %>%
  pull(Specification_ID)

if (length(primary_tgarch_x_id) == 1L) {
  primary_x_coverage <- results_df %>%
    filter(Specification_ID == primary_tgarch_x_id) %>%
    summarise(
      Valid_Coverage = mean(Forecast_Valid, na.rm = TRUE)
    ) %>%
    pull(Valid_Coverage)
  
  if (
    is.finite(primary_x_coverage) &&
    primary_x_coverage < MIN_TGARCH_X_VALID_COVERAGE
  ) {
    warning(
      sprintf(
        paste0(
          "Primary TGARCH-X valid coverage is %.2f%%, below the %.2f%% ",
          "review threshold. Inspect 13A_TGARCH_X_Imputation_Audit.csv ",
          "and 13_Failure_Reasons.csv before exporting to Python."
        ),
        100 * primary_x_coverage,
        100 * MIN_TGARCH_X_VALID_COVERAGE
      )
    )
  }
}

# =============================================================================
# 8. WEIGHTS, LOSSES, AND ACCURACY
# =============================================================================
losses_df <- results_df %>%
  left_join(weight_panel, by = c("Date", "Sector", "Sector_Cap_Lag1")) %>%
  mutate(
    Forecast_Positive = ifelse(
      is.finite(Forecast_Variance) & Forecast_Variance > 0,
      pmax(Forecast_Variance, FORECAST_EPSILON),
      NA_real_
    ),
    Operational_Forecast_Positive = ifelse(
      is.finite(Operational_Forecast_Variance) &
        Operational_Forecast_Variance > 0,
      pmax(Operational_Forecast_Variance, FORECAST_EPSILON),
      NA_real_
    ),
    Actual_Positive = pmax(Actual_Variance, TARGET_EPSILON),
    QLIKE_Ratio = Actual_Positive / Forecast_Positive,
    QLIKE = QLIKE_Ratio - log(QLIKE_Ratio) - 1,
    Squared_Error = (Actual_Positive - Forecast_Positive)^2,
    Absolute_Error = abs(Actual_Positive - Forecast_Positive),
    Operational_QLIKE_Ratio =
      Actual_Positive / Operational_Forecast_Positive,
    Operational_QLIKE =
      Operational_QLIKE_Ratio - log(Operational_QLIKE_Ratio) - 1,
    Operational_Squared_Error =
      (Actual_Positive - Operational_Forecast_Positive)^2,
    Operational_Absolute_Error =
      abs(Actual_Positive - Operational_Forecast_Positive),
    Parkinson_QLIKE_Ratio = Parkinson_Variance / Forecast_Positive,
    Parkinson_QLIKE = ifelse(
      is.finite(Parkinson_QLIKE_Ratio) & Parkinson_QLIKE_Ratio > 0,
      Parkinson_QLIKE_Ratio - log(Parkinson_QLIKE_Ratio) - 1,
      NA_real_
    )
  )

common_support <- losses_df %>%
  group_by(Sector, Date, Evaluation_Phase) %>%
  summarise(
    N_Specifications_Valid = sum(is.finite(Forecast_Positive)),
    N_Specifications_Expected = nrow(specification_grid),
    Common_Support =
      N_Specifications_Valid == N_Specifications_Expected,
    .groups = "drop"
  )

losses_df <- losses_df %>%
  left_join(
    common_support,
    by = c("Sector", "Date", "Evaluation_Phase")
  )

sector_metrics <- losses_df %>%
  filter(is.finite(QLIKE)) %>%
  group_by(
    Specification_ID,
    Model_Family,
    Model_Label,
    Analysis_Set,
    Window_Observations,
    Distribution,
    X_Transformation,
    X_Lag,
    Sector,
    Evaluation_Phase
  ) %>%
  summarise(
    N = n(),
    RMSE = sqrt(mean(Squared_Error)),
    MAE = mean(Absolute_Error),
    QLIKE = mean(QLIKE),
    Mean_Actual_Variance = mean(Actual_Positive),
    .groups = "drop"
  )

weights_long <- losses_df %>%
  pivot_longer(
    cols = c(
      Market_Weight_Daily,
      Equal_Weight_Daily,
      Capped_Market_Weight_Daily
    ),
    names_to = "Weight_Scheme",
    values_to = "Evaluation_Weight"
  )

daily_weighted_losses <- weights_long %>%
  filter(
    is.finite(QLIKE),
    is.finite(Evaluation_Weight),
    Evaluation_Weight >= 0
  ) %>%
  group_by(
    Specification_ID,
    Model_Family,
    Model_Label,
    Analysis_Set,
    Window_Observations,
    Distribution,
    X_Transformation,
    X_Lag,
    Weight_Scheme,
    Date,
    Evaluation_Phase
  ) %>%
  summarise(
    Weight_Sum = sum(Evaluation_Weight),
    Weighted_MSE = sum(Evaluation_Weight * Squared_Error) / Weight_Sum,
    Weighted_MAE = sum(Evaluation_Weight * Absolute_Error) / Weight_Sum,
    Weighted_QLIKE = sum(Evaluation_Weight * QLIKE) / Weight_Sum,
    Weighted_Parkinson_QLIKE = if (
      any(is.finite(Parkinson_QLIKE))
    ) {
      sum(
        Evaluation_Weight[is.finite(Parkinson_QLIKE)] *
          Parkinson_QLIKE[is.finite(Parkinson_QLIKE)]
      ) / sum(Evaluation_Weight[is.finite(Parkinson_QLIKE)])
    } else {
      NA_real_
    },
    N_Sectors = n(),
    .groups = "drop"
  )

market_metrics <- daily_weighted_losses %>%
  group_by(
    Specification_ID,
    Model_Family,
    Model_Label,
    Analysis_Set,
    Window_Observations,
    Distribution,
    X_Transformation,
    X_Lag,
    Weight_Scheme,
    Evaluation_Phase
  ) %>%
  summarise(
    N_Dates = n(),
    RMSE = sqrt(mean(Weighted_MSE)),
    MAE = mean(Weighted_MAE),
    QLIKE = mean(Weighted_QLIKE),
    Parkinson_QLIKE = if (
      any(is.finite(Weighted_Parkinson_QLIKE))
    ) {
      mean(Weighted_Parkinson_QLIKE, na.rm = TRUE)
    } else {
      NA_real_
    },
    Mean_N_Sectors = mean(N_Sectors),
    .groups = "drop"
  )

operational_metrics <- losses_df %>%
  filter(is.finite(Operational_QLIKE)) %>%
  group_by(
    Specification_ID,
    Model_Family,
    Sector,
    Evaluation_Phase
  ) %>%
  summarise(
    N = n(),
    N_Fallback = sum(Operational_Fallback_Used, na.rm = TRUE),
    Fallback_Rate = mean(Operational_Fallback_Used, na.rm = TRUE),
    RMSE = sqrt(mean(Operational_Squared_Error)),
    MAE = mean(Operational_Absolute_Error),
    QLIKE = mean(Operational_QLIKE),
    .groups = "drop"
  )

write_csv_safe(losses_df, "14_Forecasts_With_Losses.csv")
write_csv_safe(common_support, "15_Common_Support_Map.csv")
write_csv_safe(sector_metrics, "16_Accuracy_by_Sector.csv")
write_csv_safe(daily_weighted_losses, "17_Daily_Weighted_Losses.csv")
write_csv_safe(market_metrics, "18_Accuracy_by_Weight_Scheme.csv")
write_csv_safe(operational_metrics, "19_Operational_Accuracy_With_Fallback.csv")

# =============================================================================
# 9. PYTHON-READY EXPORTS
# =============================================================================
common_target_python <- common_actual %>%
  left_join(weight_panel, by = c("Date", "Sector", "Sector_Cap_Lag1")) %>%
  arrange(Sector, Date)

write_csv_safe(
  common_target_python,
  "20_Common_Target_for_Python.csv"
)

if (RUN_PROFILE %in% c("PRIMARY", "ALL")) {
  primary_results <- results_df %>%
    filter(Primary_Flag)
  
  primary_model_names <- c(
    "RollingVariance" = "Historical_Variance_Forecast",
    "EWMA" = "EWMA_Forecast",
    "sGARCH" = "GARCH_Forecast",
    "gjrGARCH" = "GJR_GARCH_Forecast",
    "fGARCH-TGARCH" = "TGARCH_Forecast_NoX",
    "fGARCH-TGARCH-X" = "TGARCH_X_Forecast"
  )
  
  primary_wide <- primary_results %>%
    mutate(
      Python_Column = unname(primary_model_names[Model_Family])
    ) %>%
    select(
      Date,
      Sector,
      Python_Column,
      Forecast_Variance
    ) %>%
    pivot_wider(
      names_from = Python_Column,
      values_from = Forecast_Variance
    ) %>%
    left_join(
      common_target_python,
      by = c("Date", "Sector")
    ) %>%
    arrange(Sector, Date)
  
  write_csv_safe(
    primary_wide,
    "21_Econometric_Primary_Forecasts_Wide_for_Python.csv"
  )
  
  focal_tgarch_x <- primary_results %>%
    filter(
      Model_Family == "fGARCH-TGARCH-X",
      Window_Observations == PRIMARY_WINDOW_SIZE,
      Distribution == PRIMARY_DISTRIBUTION,
      X_Transformation == PRIMARY_X_TRANSFORMATION,
      X_Lag == PRIMARY_X_LAG,
      Forecast_Valid
    ) %>%
    left_join(
      weight_panel,
      by = c("Date", "Sector", "Sector_Cap_Lag1")
    ) %>%
    transmute(
      Date,
      Sector,
      TGARCH_Forecast = Forecast_Variance,
      Actual_Variance,
      Actual_Variance_Raw,
      Sector_Return,
      Sector_Cap_Lag1,
      Market_Weight_Daily,
      Equal_Weight_Daily,
      Capped_Market_Weight_Daily,
      N_Stocks,
      Volume_Coverage,
      Range_Coverage,
      Parkinson_Variance,
      Window_Size = Window_Observations,
      Distribution,
      X_Transformation,
      X_Lag,
      Forecast_Valid,
      Operational_Fallback_Used,
      Code_Version
    ) %>%
    arrange(Sector, Date)
  
  if (nrow(focal_tgarch_x) == 0L) {
    stop("No valid primary TGARCH-X forecasts were available for Python.")
  }
  
  if (
    any(!is.finite(focal_tgarch_x$Market_Weight_Daily)) ||
    any(!is.finite(focal_tgarch_x$Equal_Weight_Daily)) ||
    any(!is.finite(focal_tgarch_x$Capped_Market_Weight_Daily))
  ) {
    stop("The focal TGARCH-X export contains missing evaluation weights.")
  }
  
  if (anyDuplicated(focal_tgarch_x[, c("Sector", "Date")]) > 0L) {
    stop("Duplicate rows exist in the clean TGARCH-X Python export.")
  }
  
  write_csv_safe(
    focal_tgarch_x,
    "TGARCH_X_Forecasts_CLEAN.csv"
  )
}

# =============================================================================
# 10. PAIRWISE QLIKE TESTS WITH HAC, BLOCK BOOTSTRAP, AND HOLM ADJUSTMENT
# =============================================================================
automatic_hac_lag <- function(n) {
  max(0L, floor(4 * (n / 100)^(2 / 9)))
}

newey_west_long_run_variance <- function(x, lag_value) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L) return(NA_real_)
  
  centered <- x - mean(x)
  long_run_variance <- sum(centered^2) / n
  
  if (lag_value > 0L) {
    for (k in seq_len(min(lag_value, n - 1L))) {
      covariance_k <- sum(
        centered[(k + 1L):n] * centered[1L:(n - k)]
      ) / n
      bartlett_weight <- 1 - k / (lag_value + 1)
      long_run_variance <- long_run_variance +
        2 * bartlett_weight * covariance_k
    }
  }
  long_run_variance
}

moving_block_bootstrap_ci <- function(
    x,
    replications = BOOTSTRAP_REPLICATIONS,
    block_length = BOOTSTRAP_BLOCK_LENGTH,
    confidence = 0.95
) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 20L) {
    return(data.frame(
      Bootstrap_Block_Length = NA_integer_,
      Bootstrap_CI_Lower = NA_real_,
      Bootstrap_CI_Upper = NA_real_
    ))
  }
  
  block_value <- if (is.na(block_length)) {
    max(2L, ceiling(n^(1 / 3)))
  } else {
    max(2L, min(as.integer(block_length), n))
  }
  
  number_of_blocks <- ceiling(n / block_value)
  bootstrap_means <- numeric(replications)
  
  for (b in seq_len(replications)) {
    starts <- sample.int(n, number_of_blocks, replace = TRUE)
    indices <- unlist(
      lapply(starts, function(start) {
        ((start - 1L + seq_len(block_value) - 1L) %% n) + 1L
      }),
      use.names = FALSE
    )
    indices <- indices[seq_len(n)]
    bootstrap_means[[b]] <- mean(x[indices])
  }
  
  alpha <- 1 - confidence
  quantiles <- quantile(
    bootstrap_means,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE,
    type = 6
  )
  
  data.frame(
    Bootstrap_Block_Length = block_value,
    Bootstrap_CI_Lower = quantiles[[1L]],
    Bootstrap_CI_Upper = quantiles[[2L]]
  )
}

dm_hac_test <- function(loss_a, loss_b, horizon = 1L, hac_lag = NA_integer_) {
  complete <- is.finite(loss_a) & is.finite(loss_b)
  differential <- loss_a[complete] - loss_b[complete]
  n <- length(differential)
  
  if (n < 20L) {
    return(data.frame(
      N = n,
      Mean_Loss_Difference = if (n > 0L) mean(differential) else NA_real_,
      HAC_Lag = NA_integer_,
      Statistic = NA_real_,
      P_Value = NA_real_,
      Bootstrap_Block_Length = NA_integer_,
      Bootstrap_CI_Lower = NA_real_,
      Bootstrap_CI_Upper = NA_real_
    ))
  }
  
  lag_value <- if (is.na(hac_lag)) automatic_hac_lag(n) else max(0L, as.integer(hac_lag))
  lrv <- newey_west_long_run_variance(differential, lag_value)
  bootstrap_ci <- moving_block_bootstrap_ci(differential)
  
  if (!is.finite(lrv) || lrv <= 0) {
    return(cbind(
      data.frame(
        N = n,
        Mean_Loss_Difference = mean(differential),
        HAC_Lag = lag_value,
        Statistic = NA_real_,
        P_Value = NA_real_
      ),
      bootstrap_ci
    ))
  }
  
  statistic <- mean(differential) / sqrt(lrv / n)
  hln_factor <- sqrt(
    (n + 1 - 2 * horizon + horizon * (horizon - 1) / n) / n
  )
  statistic_hln <- statistic * hln_factor
  p_value <- 2 * pt(abs(statistic_hln), df = n - 1L, lower.tail = FALSE)
  
  cbind(
    data.frame(
      N = n,
      Mean_Loss_Difference = mean(differential),
      HAC_Lag = lag_value,
      Statistic = statistic_hln,
      P_Value = p_value
    ),
    bootstrap_ci
  )
}

pairwise_tests <- data.frame(
  Message = "Pairwise primary tests are produced only in PRIMARY or ALL profiles."
)

if (RUN_PROFILE %in% c("PRIMARY", "ALL")) {
  primary_losses <- losses_df %>%
    filter(Primary_Flag, is.finite(QLIKE)) %>%
    select(Sector, Date, Evaluation_Phase, Model_Family, QLIKE) %>%
    pivot_wider(names_from = Model_Family, values_from = QLIKE)
  
  target_model <- "fGARCH-TGARCH-X"
  comparators <- setdiff(
    c(
      "RollingVariance",
      "EWMA",
      "sGARCH",
      "gjrGARCH",
      "fGARCH-TGARCH"
    ),
    target_model
  )
  
  pairwise_rows <- list()
  pairwise_index <- 0L
  set.seed(SEED)
  
  for (phase_name in unique(primary_losses$Evaluation_Phase)) {
    phase_data <- primary_losses %>%
      filter(Evaluation_Phase == phase_name)
    
    for (sector_name in unique(phase_data$Sector)) {
      sector_data <- phase_data %>%
        filter(Sector == sector_name) %>%
        arrange(Date)
      
      for (comparator in comparators) {
        if (
          !target_model %in% names(sector_data) ||
          !comparator %in% names(sector_data)
        ) next
        
        test_result <- dm_hac_test(
          loss_a = sector_data[[target_model]],
          loss_b = sector_data[[comparator]],
          horizon = FORECAST_HORIZON,
          hac_lag = DM_HAC_LAG
        )
        
        pairwise_index <- pairwise_index + 1L
        pairwise_rows[[pairwise_index]] <- data.frame(
          Sector = sector_name,
          Evaluation_Phase = phase_name,
          Target_Model = target_model,
          Comparator = comparator,
          test_result,
          Direction = ifelse(
            is.finite(test_result$Mean_Loss_Difference) &
              test_result$Mean_Loss_Difference < 0,
            "TGARCH-X lower QLIKE",
            "Comparator lower or equal QLIKE"
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  pairwise_tests <- if (length(pairwise_rows) > 0L) {
    bind_rows(pairwise_rows) %>%
      group_by(Evaluation_Phase) %>%
      mutate(
        P_Value_Holm = p.adjust(P_Value, method = P_VALUE_ADJUSTMENT),
        P_Value_BH = p.adjust(P_Value, method = "BH")
      ) %>%
      ungroup()
  } else {
    data.frame(Message = "No primary pairwise tests could be computed.")
  }
}

write_csv_safe(pairwise_tests, "22_Pairwise_QLIKE_HAC_Bootstrap_Tests.csv")

# =============================================================================
# 11. RESIDUAL DIAGNOSTICS FOR THE LAST SUCCESSFUL FIT OF EACH SPECIFICATION
# =============================================================================
arch_lm_test <- function(standardized_residuals, lags = 10L) {
  z <- as.numeric(standardized_residuals)
  z <- z[is.finite(z)]
  if (length(z) <= lags + 5L) {
    return(data.frame(Statistic = NA_real_, DF = lags, P_Value = NA_real_))
  }
  
  squared <- z^2
  embedded <- embed(squared, lags + 1L)
  dependent <- embedded[, 1L]
  regressors <- embedded[, -1L, drop = FALSE]
  auxiliary_fit <- lm(dependent ~ regressors)
  statistic <- nrow(embedded) * summary(auxiliary_fit)$r.squared
  
  data.frame(
    Statistic = statistic,
    DF = lags,
    P_Value = pchisq(statistic, df = lags, lower.tail = FALSE)
  )
}

diagnostic_rows <- list()
sign_bias_rows <- list()
diagnostic_index <- 0L
sign_bias_index <- 0L

for (fit_key in names(last_successful_fit)) {
  fit_object <- last_successful_fit[[fit_key]]
  metadata <- last_successful_metadata[[fit_key]]
  
  standardized_residuals <- tryCatch(
    as.numeric(residuals(fit_object, standardize = TRUE)),
    error = function(e) numeric(0)
  )
  standardized_residuals <- standardized_residuals[
    is.finite(standardized_residuals)
  ]
  
  diagnostic_lag <- min(
    DIAGNOSTIC_LAG,
    max(1L, floor(length(standardized_residuals) / 5L))
  )
  
  lb_residual <- if (length(standardized_residuals) > diagnostic_lag + 2L) {
    Box.test(
      standardized_residuals,
      lag = diagnostic_lag,
      type = "Ljung-Box"
    )
  } else {
    NULL
  }
  
  lb_squared <- if (length(standardized_residuals) > diagnostic_lag + 2L) {
    Box.test(
      standardized_residuals^2,
      lag = diagnostic_lag,
      type = "Ljung-Box"
    )
  } else {
    NULL
  }
  
  arch_lm <- arch_lm_test(standardized_residuals, diagnostic_lag)
  
  diagnostic_index <- diagnostic_index + 1L
  diagnostic_rows[[diagnostic_index]] <- data.frame(
    Sector = metadata$Sector,
    Specification_ID = metadata$Specification_ID,
    Model_Family = metadata$Model_Family,
    Window_Start_Date = metadata$Window_Start_Date,
    Window_End_Date = metadata$Window_End_Date,
    Forecast_Date = metadata$Forecast_Date,
    Diagnostic_Lag = diagnostic_lag,
    Ljung_Box_Residual_Statistic = if (is.null(lb_residual)) {
      NA_real_
    } else {
      as.numeric(lb_residual$statistic)
    },
    Ljung_Box_Residual_P_Value = if (is.null(lb_residual)) {
      NA_real_
    } else {
      lb_residual$p.value
    },
    Ljung_Box_Squared_Residual_Statistic = if (is.null(lb_squared)) {
      NA_real_
    } else {
      as.numeric(lb_squared$statistic)
    },
    Ljung_Box_Squared_Residual_P_Value = if (is.null(lb_squared)) {
      NA_real_
    } else {
      lb_squared$p.value
    },
    ARCH_LM_Statistic = arch_lm$Statistic,
    ARCH_LM_DF = arch_lm$DF,
    ARCH_LM_P_Value = arch_lm$P_Value,
    Persistence = tryCatch(
      as.numeric(rugarch::persistence(fit_object)),
      error = function(e) NA_real_
    ),
    Code_Version = CODE_VERSION,
    stringsAsFactors = FALSE
  )
  
  sign_bias_output <- tryCatch(
    rugarch::signbias(fit_object),
    error = function(e) NULL
  )
  if (!is.null(sign_bias_output)) {
    sign_bias_index <- sign_bias_index + 1L
    sign_bias_part <- as.data.frame(sign_bias_output, stringsAsFactors = FALSE)
    sign_bias_part$Test <- rownames(sign_bias_part)
    rownames(sign_bias_part) <- NULL
    sign_bias_part$Sector <- metadata$Sector
    sign_bias_part$Specification_ID <- metadata$Specification_ID
    sign_bias_part$Model_Family <- metadata$Model_Family
    sign_bias_part$Window_End_Date <- metadata$Window_End_Date
    sign_bias_part$Forecast_Date <- metadata$Forecast_Date
    sign_bias_part$Code_Version <- CODE_VERSION
    sign_bias_rows[[sign_bias_index]] <- sign_bias_part
  }
}

diagnostics_df <- if (length(diagnostic_rows) > 0L) {
  bind_rows(diagnostic_rows)
} else {
  data.frame(Message = "No successful GARCH fit was available for diagnostics.")
}
sign_bias_df <- if (length(sign_bias_rows) > 0L) {
  bind_rows(sign_bias_rows)
} else {
  data.frame(Message = "No successful GARCH fit was available for sign-bias tests.")
}

write_csv_safe(diagnostics_df, "23_Last_Fit_Residual_Diagnostics.csv")
write_csv_safe(sign_bias_df, "24_Last_Fit_Sign_Bias_Tests.csv")

# =============================================================================
# 12. EXPLORATORY STABILITY OF TGARCH-X RELATIVE TO TGARCH
# =============================================================================
if (
  isTRUE(RUN_STABILITY_ANALYSIS) &&
  RUN_PROFILE %in% c("PRIMARY", "ALL")
) {
  stability_source <- losses_df %>%
    filter(
      Primary_Flag,
      Model_Family %in% c("fGARCH-TGARCH", "fGARCH-TGARCH-X"),
      is.finite(QLIKE)
    ) %>%
    select(Sector, Date, Evaluation_Phase, Model_Family, QLIKE) %>%
    pivot_wider(names_from = Model_Family, values_from = QLIKE) %>%
    mutate(
      Loss_Difference_X_minus_NoX =
        .data[["fGARCH-TGARCH-X"]] -
        .data[["fGARCH-TGARCH"]]
    ) %>%
    group_by(Sector) %>%
    arrange(Date, .by_group = TRUE) %>%
    mutate(
      Rolling_Mean_Loss_Difference = zoo::rollapplyr(
        Loss_Difference_X_minus_NoX,
        width = STABILITY_ROLLING_WINDOW,
        FUN = mean,
        fill = NA_real_,
        partial = FALSE,
        na.rm = TRUE
      )
    ) %>%
    ungroup()
  
  write_csv_safe(
    stability_source,
    "25_Rolling_Relative_QLIKE_TGARCH_X_vs_TGARCH.csv"
  )
  
  breakpoint_rows <- list()
  breakpoint_index <- 0L
  
  for (sector_name in unique(stability_source$Sector)) {
    sector_stability <- stability_source %>%
      filter(
        Sector == sector_name,
        is.finite(Loss_Difference_X_minus_NoX)
      ) %>%
      arrange(Date)
    
    if (nrow(sector_stability) < 80L) next
    
    breakpoint_fit <- tryCatch(
      breakpoints(
        Loss_Difference_X_minus_NoX ~ 1,
        data = sector_stability,
        h = 0.15
      ),
      error = function(e) NULL
    )
    if (is.null(breakpoint_fit)) next
    
    bic_values <- BIC(breakpoint_fit)
    selected_breaks <- which.min(bic_values) - 1L
    
    if (selected_breaks <= 0L) {
      breakpoint_index <- breakpoint_index + 1L
      breakpoint_rows[[breakpoint_index]] <- data.frame(
        Sector = sector_name,
        Number_of_Breaks = 0L,
        Breakpoint_Index = NA_integer_,
        Breakpoint_Date = as.Date(NA),
        Code_Version = CODE_VERSION,
        stringsAsFactors = FALSE
      )
      next
    }
    
    selected_fit <- breakpoints(breakpoint_fit, breaks = selected_breaks)
    indices <- selected_fit$breakpoints
    indices <- indices[is.finite(indices)]
    
    for (break_index in indices) {
      breakpoint_index <- breakpoint_index + 1L
      breakpoint_rows[[breakpoint_index]] <- data.frame(
        Sector = sector_name,
        Number_of_Breaks = selected_breaks,
        Breakpoint_Index = break_index,
        Breakpoint_Date = sector_stability$Date[[break_index]],
        Code_Version = CODE_VERSION,
        stringsAsFactors = FALSE
      )
    }
  }
  
  breakpoint_df <- if (length(breakpoint_rows) > 0L) {
    bind_rows(breakpoint_rows)
  } else {
    data.frame(Message = "No eligible stability breakpoints were estimated.")
  }
  write_csv_safe(
    breakpoint_df,
    "26_Exploratory_Relative_Loss_Breakpoints.csv"
  )
}


# =============================================================================
# 12A. REVIEWER-ALIGNED SENSITIVITY SUMMARY ON COMMON FINAL-TEST SUPPORT
# =============================================================================
# The frozen primary R outputs remain authoritative. This section compares
# sensitivity alternatives with the re-estimated within-run primary anchor on
# exactly the same FinalTest sector-days. This avoids differences in forecast
# availability mechanically driving robustness conclusions.

is_primary_anchor <- function(window_size, distribution, x_transformation, x_lag) {
  window_size == PRIMARY_WINDOW_SIZE &
    distribution == PRIMARY_DISTRIBUTION &
    x_transformation == PRIMARY_X_TRANSFORMATION &
    x_lag == PRIMARY_X_LAG
}

specification_dimension_map <- list()
map_index <- 0L

for (i in seq_len(nrow(specification_grid))) {
  s <- specification_grid[i, , drop = FALSE]
  is_anchor <- is_primary_anchor(
    s$Window_Size[[1L]],
    s$Distribution[[1L]],
    s$X_Transformation[[1L]],
    s$X_Lag[[1L]]
  )
  
  add_map_row <- function(family, label, anchor_flag) {
    map_index <<- map_index + 1L
    specification_dimension_map[[map_index]] <<- data.frame(
      Specification_ID = s$Specification_ID[[1L]],
      Sensitivity_Family = family,
      Sensitivity_Label = label,
      Primary_Anchor = anchor_flag,
      stringsAsFactors = FALSE
    )
  }
  
  if (is_anchor) {
    add_map_row(
      "WINDOW",
      paste0(PRIMARY_WINDOW_SIZE, "-observation primary anchor"),
      TRUE
    )
    add_map_row(
      "DISTRIBUTION",
      "Student-t primary anchor",
      TRUE
    )
    add_map_row(
      "X",
      "Turnover change lag 1 primary anchor",
      TRUE
    )
  } else if (
    s$Window_Size[[1L]] != PRIMARY_WINDOW_SIZE &&
    s$Distribution[[1L]] == PRIMARY_DISTRIBUTION &&
    s$X_Transformation[[1L]] == PRIMARY_X_TRANSFORMATION &&
    s$X_Lag[[1L]] == PRIMARY_X_LAG
  ) {
    add_map_row(
      "WINDOW",
      paste0(s$Window_Size[[1L]], "-observation window"),
      FALSE
    )
  } else if (
    s$Window_Size[[1L]] == PRIMARY_WINDOW_SIZE &&
    s$Distribution[[1L]] != PRIMARY_DISTRIBUTION &&
    s$X_Transformation[[1L]] == PRIMARY_X_TRANSFORMATION &&
    s$X_Lag[[1L]] == PRIMARY_X_LAG
  ) {
    distribution_label <- switch(
      s$Distribution[[1L]],
      "norm" = "Gaussian innovations",
      "std" = "Student-t innovations",
      "ged" = "GED innovations",
      paste0("Innovations: ", s$Distribution[[1L]])
    )
    add_map_row("DISTRIBUTION", distribution_label, FALSE)
  } else if (
    s$Window_Size[[1L]] == PRIMARY_WINDOW_SIZE &&
    s$Distribution[[1L]] == PRIMARY_DISTRIBUTION &&
    (
      s$X_Transformation[[1L]] != PRIMARY_X_TRANSFORMATION ||
      s$X_Lag[[1L]] != PRIMARY_X_LAG
    )
  ) {
    x_label <- if (
      s$X_Transformation[[1L]] == "turnover_change" &&
      s$X_Lag[[1L]] == 2L
    ) {
      "Turnover change lag 2"
    } else if (
      s$X_Transformation[[1L]] == "sector_volume_change" &&
      s$X_Lag[[1L]] == 1L
    ) {
      "Sector volume change lag 1"
    } else {
      paste0(
        s$X_Transformation[[1L]],
        " lag ",
        s$X_Lag[[1L]]
      )
    }
    add_map_row("X", x_label, FALSE)
  }
}

specification_dimension_map <- bind_rows(specification_dimension_map)
write_csv_safe(
  specification_dimension_map,
  "34_Sensitivity_Specification_Map.csv"
)

final_common_losses <- losses_df %>%
  filter(
    Date >= SENSITIVITY_FINAL_TEST_START,
    Date <= SAMPLE_END_DATE,
    Common_Support,
    is.finite(QLIKE)
  )

if (nrow(final_common_losses) == 0L) {
  stop(
    "No common-support FinalTest sensitivity observations are available. ",
    "Inspect 12_Fit_Coverage_and_Quality.csv and 15_Common_Support_Map.csv."
  )
}

weight_columns <- c(
  "Market_Weight_Daily",
  "Equal_Weight_Daily",
  "Capped_Market_Weight_Daily"
)

sensitivity_weights_long <- final_common_losses %>%
  pivot_longer(
    cols = all_of(weight_columns),
    names_to = "Weight_Scheme",
    values_to = "Evaluation_Weight"
  ) %>%
  filter(
    is.finite(Evaluation_Weight),
    Evaluation_Weight >= 0
  )

sensitivity_daily_losses <- sensitivity_weights_long %>%
  group_by(
    Specification_ID,
    Weight_Scheme,
    Date
  ) %>%
  summarise(
    Weight_Sum = sum(Evaluation_Weight),
    Weighted_QLIKE = sum(Evaluation_Weight * QLIKE) / Weight_Sum,
    Weighted_MSE = sum(Evaluation_Weight * Squared_Error) / Weight_Sum,
    Weighted_MAE = sum(Evaluation_Weight * Absolute_Error) / Weight_Sum,
    Weighted_Parkinson_QLIKE = if (
      any(is.finite(Parkinson_QLIKE))
    ) {
      sum(
        Evaluation_Weight[is.finite(Parkinson_QLIKE)] *
          Parkinson_QLIKE[is.finite(Parkinson_QLIKE)]
      ) / sum(Evaluation_Weight[is.finite(Parkinson_QLIKE)])
    } else {
      NA_real_
    },
    N_Sectors = n(),
    .groups = "drop"
  )

sensitivity_accuracy_core <- sensitivity_daily_losses %>%
  group_by(
    Specification_ID,
    Weight_Scheme
  ) %>%
  summarise(
    N_Dates = n(),
    Mean_N_Sectors = mean(N_Sectors),
    QLIKE = mean(Weighted_QLIKE),
    RMSE = sqrt(mean(Weighted_MSE)),
    MAE = mean(Weighted_MAE),
    Parkinson_QLIKE = if (
      any(is.finite(Weighted_Parkinson_QLIKE))
    ) {
      mean(Weighted_Parkinson_QLIKE, na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

sensitivity_accuracy <- specification_dimension_map %>%
  left_join(
    sensitivity_accuracy_core,
    by = "Specification_ID"
  ) %>%
  group_by(
    Sensitivity_Family,
    Weight_Scheme
  ) %>%
  mutate(
    Anchor_QLIKE = QLIKE[Primary_Anchor][1L],
    Delta_QLIKE_vs_Anchor = QLIKE - Anchor_QLIKE,
    Anchor_RMSE = RMSE[Primary_Anchor][1L],
    Delta_RMSE_vs_Anchor = RMSE - Anchor_RMSE,
    Anchor_MAE = MAE[Primary_Anchor][1L],
    Delta_MAE_vs_Anchor = MAE - Anchor_MAE,
    QLIKE_Rank = rank(QLIKE, ties.method = "min"),
    Direction_vs_Anchor = case_when(
      Primary_Anchor ~ "Primary anchor",
      Delta_QLIKE_vs_Anchor < 0 ~ "Alternative lower QLIKE",
      Delta_QLIKE_vs_Anchor > 0 ~ "Primary anchor lower QLIKE",
      TRUE ~ "Equal QLIKE"
    )
  ) %>%
  ungroup() %>%
  left_join(
    specification_grid %>%
      select(
        Specification_ID,
        Window_Size,
        Distribution,
        X_Transformation,
        X_Lag
      ),
    by = "Specification_ID"
  ) %>%
  arrange(
    Sensitivity_Family,
    Weight_Scheme,
    QLIKE_Rank,
    Sensitivity_Label
  )

write_csv_safe(
  sensitivity_accuracy,
  "35_Sensitivity_CommonSupport_FinalTest_Accuracy.csv"
)

# Sector-level QLIKE on the same common support.
sensitivity_sector_accuracy <- final_common_losses %>%
  group_by(
    Specification_ID,
    Sector
  ) %>%
  summarise(
    N = n(),
    QLIKE = mean(QLIKE),
    RMSE = sqrt(mean(Squared_Error)),
    MAE = mean(Absolute_Error),
    .groups = "drop"
  ) %>%
  inner_join(
    specification_dimension_map,
    by = "Specification_ID"
  ) %>%
  group_by(
    Sensitivity_Family,
    Sector
  ) %>%
  mutate(
    Anchor_QLIKE = QLIKE[Primary_Anchor][1L],
    Delta_QLIKE_vs_Anchor = QLIKE - Anchor_QLIKE
  ) %>%
  ungroup() %>%
  arrange(
    Sensitivity_Family,
    Sector,
    Primary_Anchor,
    Sensitivity_Label
  )

write_csv_safe(
  sensitivity_sector_accuracy,
  "36_Sensitivity_CommonSupport_FinalTest_Sector_Accuracy.csv"
)

# -------------------------------------------------------------------------
# Paired daily market-capitalisation QLIKE tests against the anchor
# -------------------------------------------------------------------------
market_daily <- sensitivity_daily_losses %>%
  filter(Weight_Scheme == "Market_Weight_Daily") %>%
  select(
    Date,
    Specification_ID,
    Weighted_QLIKE
  )

anchor_id <- specification_grid %>%
  filter(
    Model_Family == "fGARCH-TGARCH-X",
    Window_Size == PRIMARY_WINDOW_SIZE,
    Distribution == PRIMARY_DISTRIBUTION,
    X_Transformation == PRIMARY_X_TRANSFORMATION,
    X_Lag == PRIMARY_X_LAG
  ) %>%
  pull(Specification_ID)

if (length(anchor_id) != 1L) {
  stop("Sensitivity suite could not identify exactly one primary-anchor specification.")
}

anchor_daily <- market_daily %>%
  filter(Specification_ID == anchor_id[[1L]]) %>%
  select(
    Date,
    Anchor_QLIKE = Weighted_QLIKE
  )

pairwise_sensitivity_rows <- list()
pairwise_sensitivity_index <- 0L
set.seed(SEED)

alternative_map <- specification_dimension_map %>%
  filter(!Primary_Anchor)

for (i in seq_len(nrow(alternative_map))) {
  alt <- alternative_map[i, , drop = FALSE]
  paired <- market_daily %>%
    filter(Specification_ID == alt$Specification_ID[[1L]]) %>%
    select(
      Date,
      Alternative_QLIKE = Weighted_QLIKE
    ) %>%
    inner_join(anchor_daily, by = "Date") %>%
    arrange(Date)
  
  test_result <- dm_hac_test(
    loss_a = paired$Alternative_QLIKE,
    loss_b = paired$Anchor_QLIKE,
    horizon = FORECAST_HORIZON,
    hac_lag = DM_HAC_LAG
  )
  
  pairwise_sensitivity_index <- pairwise_sensitivity_index + 1L
  pairwise_sensitivity_rows[[pairwise_sensitivity_index]] <- data.frame(
    Sensitivity_Family = alt$Sensitivity_Family[[1L]],
    Alternative_Label = alt$Sensitivity_Label[[1L]],
    Alternative_Specification_ID = alt$Specification_ID[[1L]],
    Anchor_Specification_ID = anchor_id[[1L]],
    Weight_Scheme = "Market_Weight_Daily",
    Date_Start = if (nrow(paired) > 0L) min(paired$Date) else as.Date(NA),
    Date_End = if (nrow(paired) > 0L) max(paired$Date) else as.Date(NA),
    test_result,
    Sign_Convention = "Negative mean difference favours alternative",
    stringsAsFactors = FALSE
  )
}

sensitivity_pairwise_tests <- bind_rows(pairwise_sensitivity_rows) %>%
  group_by(Sensitivity_Family) %>%
  mutate(
    P_Value_Holm_Within_Family = p.adjust(P_Value, method = "holm")
  ) %>%
  ungroup() %>%
  mutate(
    P_Value_Holm_Global = p.adjust(P_Value, method = "holm"),
    P_Value_BH_Global = p.adjust(P_Value, method = "BH"),
    Direction = case_when(
      is.finite(Mean_Loss_Difference) & Mean_Loss_Difference < 0 ~
        "Alternative lower QLIKE",
      is.finite(Mean_Loss_Difference) & Mean_Loss_Difference > 0 ~
        "Primary anchor lower QLIKE",
      TRUE ~ "Equal/undefined"
    )
  ) %>%
  arrange(
    Sensitivity_Family,
    Alternative_Label
  )

write_csv_safe(
  sensitivity_pairwise_tests,
  "37_Sensitivity_Paired_QLIKE_HAC_HLN_Bootstrap_Holm.csv"
)

# -------------------------------------------------------------------------
# Fit coverage and common-support audit
# -------------------------------------------------------------------------
sensitivity_coverage <- results_df %>%
  group_by(
    Specification_ID,
    Window_Observations,
    Distribution,
    X_Transformation,
    X_Lag
  ) %>%
  summarise(
    N_Origins = n(),
    N_Valid = sum(Forecast_Valid, na.rm = TRUE),
    Valid_Coverage = mean(Forecast_Valid, na.rm = TRUE),
    N_Operational_Fallback = sum(Operational_Fallback_Used, na.rm = TRUE),
    Operational_Fallback_Rate = mean(
      Operational_Fallback_Used,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    specification_dimension_map,
    by = "Specification_ID"
  ) %>%
  arrange(
    Sensitivity_Family,
    Primary_Anchor,
    Sensitivity_Label
  )

write_csv_safe(
  sensitivity_coverage,
  "38_Sensitivity_Fit_Coverage_Summary.csv"
)

common_support_audit <- common_support %>%
  mutate(
    FinalTest = Date >= SENSITIVITY_FINAL_TEST_START &
      Date <= SAMPLE_END_DATE
  ) %>%
  group_by(FinalTest) %>%
  summarise(
    N_Sector_Days = n(),
    N_Common_Support = sum(Common_Support, na.rm = TRUE),
    Common_Support_Rate = mean(Common_Support, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_safe(
  common_support_audit,
  "39_Sensitivity_Common_Support_Audit.csv"
)

# Compact reviewer-facing conclusion table. These labels describe direction
# only; statistical significance is taken from file 37.
market_cap_summary <- sensitivity_accuracy %>%
  filter(Weight_Scheme == "Market_Weight_Daily") %>%
  select(
    Sensitivity_Family,
    Sensitivity_Label,
    Primary_Anchor,
    QLIKE,
    Delta_QLIKE_vs_Anchor,
    QLIKE_Rank,
    Parkinson_QLIKE,
    N_Dates,
    Mean_N_Sectors
  ) %>%
  left_join(
    sensitivity_pairwise_tests %>%
      select(
        Sensitivity_Family,
        Alternative_Label,
        Mean_Loss_Difference,
        P_Value,
        P_Value_Holm_Within_Family,
        P_Value_Holm_Global,
        Bootstrap_CI_Lower,
        Bootstrap_CI_Upper
      ),
    by = c(
      "Sensitivity_Family",
      "Sensitivity_Label" = "Alternative_Label"
    )
  ) %>%
  mutate(
    Reviewer_Interpretation = case_when(
      Primary_Anchor ~ "Primary anchor",
      Delta_QLIKE_vs_Anchor < 0 &
        is.finite(P_Value_Holm_Global) &
        P_Value_Holm_Global < 0.05 ~
        "Alternative significantly improves QLIKE after global Holm correction",
      Delta_QLIKE_vs_Anchor < 0 ~
        "Alternative lowers QLIKE but not robust after global Holm correction",
      Delta_QLIKE_vs_Anchor > 0 &
        is.finite(P_Value_Holm_Global) &
        P_Value_Holm_Global < 0.05 ~
        "Alternative significantly worsens QLIKE after global Holm correction",
      Delta_QLIKE_vs_Anchor > 0 ~
        "Alternative raises QLIKE; difference not robust after global Holm correction",
      TRUE ~ "No material directional difference"
    )
  ) %>%
  arrange(
    Sensitivity_Family,
    Primary_Anchor,
    Sensitivity_Label
  )

write_csv_safe(
  market_cap_summary,
  "40_Sensitivity_Reviewer_Conclusion_Summary.csv"
)

# =============================================================================
# 13. AUDITS, PARAMETERIZATION NOTE, HASHES, AND FINAL SUMMARY
# =============================================================================
study_period_audit <- data.frame(
  Requested_Start = SAMPLE_START_DATE,
  Requested_End = SAMPLE_END_DATE,
  Observed_Start_After_Cleaning = min(df$Date),
  Observed_End_After_Cleaning = max(df$Date),
  N_Firms = n_distinct(df$Stock),
  N_Sectors = n_distinct(df$Sector),
  stringsAsFactors = FALSE
)
write_csv_safe(study_period_audit, "27_Study_Period_Audit.csv")

sector_coverage_audit <- df %>%
  group_by(Sector) %>%
  summarise(
    N_Firms = n_distinct(Stock),
    N_Rows = n(),
    N_Dates = n_distinct(Date),
    Observed_Start = min(Date),
    Observed_End = max(Date),
    High_Low_Coverage = mean(
      is.finite(High) & is.finite(Low) & High > 0 & Low > 0 & High >= Low
    ),
    Volume_Coverage = mean(is.finite(Trading_Volume) & Trading_Volume >= 0),
    .groups = "drop"
  ) %>%
  arrange(Sector)

write_csv_safe(sector_coverage_audit, "27A_Sector_Coverage_Audit.csv")

if (n_distinct(df$Stock) != EXPECTED_N_FIRMS) {
  warning(
    "Expected ", EXPECTED_N_FIRMS, " firms after core cleaning, but found ",
    n_distinct(df$Stock), ". Check the firm-inclusion audit."
  )
}
if (n_distinct(df$Sector) != EXPECTED_N_SECTORS) {
  warning(
    "Expected ", EXPECTED_N_SECTORS, " sectors after core cleaning, but found ",
    n_distinct(df$Sector), ". Check sector labels and classifications."
  )
}

if (min(df$Date) > SAMPLE_START_DATE + 31L) {
  warning(
    "The observed cleaned data begin materially after SAMPLE_START_DATE. ",
    "Do not report the requested start date until the input file is reconciled."
  )
}
if (max(df$Date) < SAMPLE_END_DATE - 31L) {
  warning(
    "The observed cleaned data end materially before SAMPLE_END_DATE. ",
    "Do not report the requested end date until the input file is reconciled."
  )
}

observation_flow <- data.frame(
  Stage = c(
    "Raw rows",
    "Invalid core rows",
    "Valid firm-price rows in study period",
    "Valid one-trading-day firm returns",
    "Sector-day target rows",
    "Rolling forecast rows",
    "Valid rolling forecasts",
    "Valid focal TGARCH-X rows for Python"
  ),
  N = c(
    nrow(raw_df),
    nrow(invalid_core_rows),
    nrow(df),
    nrow(firm_panel),
    nrow(common_actual),
    nrow(results_df),
    sum(results_df$Forecast_Valid, na.rm = TRUE),
    if (exists("focal_tgarch_x")) nrow(focal_tgarch_x) else NA_integer_
  ),
  stringsAsFactors = FALSE
)
write_csv_safe(observation_flow, "28_Observation_Flow.csv")

configuration <- data.frame(
  Setting = c(
    "CODE_VERSION",
    "PARENT_ENGINE_VERSION",
    "PARENT_ENGINE_SHA256",
    "RUN_PROFILE",
    "SEED",
    "SAMPLE_START_DATE",
    "SAMPLE_END_DATE",
    "EXPECTED_N_FIRMS",
    "EXPECTED_N_SECTORS",
    "PERIOD_TAG",
    "PRIMARY_WINDOW_SIZE",
    "WINDOW_SENSITIVITY",
    "PRIMARY_DISTRIBUTION",
    "DISTRIBUTION_SENSITIVITY",
    "RETURN_SCALE",
    "FORECAST_HORIZON",
    "EWMA_LAMBDA",
    "TARGET_EPSILON",
    "FORECAST_EPSILON",
    "MIN_VOLUME_COVERAGE",
    "MIN_RANGE_COVERAGE",
    "CAPPED_WEIGHT_MAX",
    "REQUIRE_CONTIGUOUS_SECTOR_CALENDAR",
    "CALENDAR_GAP_JUSTIFICATION",
    "PRIMARY_X_TRANSFORMATION",
    "PRIMARY_X_LAG",
    "ALLOW_X_IMPUTATION",
    "DM_HAC_LAG",
    "P_VALUE_ADJUSTMENT",
    "BOOTSTRAP_REPLICATIONS",
    "BOOTSTRAP_BLOCK_LENGTH",
    "MIN_TGARCH_X_VALID_COVERAGE",
    "RUN_STABILITY_ANALYSIS",
    "STABILITY_ROLLING_WINDOW",
    "SENSITIVITY_FINAL_TEST_START",
    "SENSITIVITY_FAMILIES",
    "MASTER_FILE"
  ),
  Value = c(
    CODE_VERSION,
    PARENT_ENGINE_VERSION,
    PARENT_ENGINE_SHA256,
    RUN_PROFILE,
    SEED,
    as.character(SAMPLE_START_DATE),
    as.character(SAMPLE_END_DATE),
    EXPECTED_N_FIRMS,
    EXPECTED_N_SECTORS,
    PERIOD_TAG,
    PRIMARY_WINDOW_SIZE,
    paste(WINDOW_SENSITIVITY, collapse = ","),
    PRIMARY_DISTRIBUTION,
    paste(DISTRIBUTION_SENSITIVITY, collapse = ","),
    RETURN_SCALE,
    FORECAST_HORIZON,
    EWMA_LAMBDA,
    TARGET_EPSILON,
    FORECAST_EPSILON,
    MIN_VOLUME_COVERAGE,
    MIN_RANGE_COVERAGE,
    CAPPED_WEIGHT_MAX,
    REQUIRE_CONTIGUOUS_SECTOR_CALENDAR,
    CALENDAR_GAP_JUSTIFICATION,
    PRIMARY_X_TRANSFORMATION,
    PRIMARY_X_LAG,
    ALLOW_X_IMPUTATION,
    DM_HAC_LAG,
    P_VALUE_ADJUSTMENT,
    BOOTSTRAP_REPLICATIONS,
    BOOTSTRAP_BLOCK_LENGTH,
    MIN_TGARCH_X_VALID_COVERAGE,
    RUN_STABILITY_ANALYSIS,
    STABILITY_ROLLING_WINDOW,
    as.character(SENSITIVITY_FINAL_TEST_START),
    paste(SENSITIVITY_FAMILIES, collapse = ","),
    normalizePath(master_file, winslash = "/", mustWork = TRUE)
  ),
  stringsAsFactors = FALSE
)
write_csv_safe(configuration, "29_Run_Configuration.csv")

parameterization_note <- c(
  "Implemented focal model: rugarch fGARCH submodel='TGARCH' with an external variance regressor.",
  "The manuscript must describe the conditional-scale parameterization implemented by rugarch.",
  "For one lag and one variance regressor, the scale recursion is represented as:",
  "sigma_t = omega + alpha1*(abs(epsilon_{t-1}) - eta11*epsilon_{t-1}) + beta1*sigma_{t-1} + vxreg1*X_{t-1}.",
  "The code standardizes X using each rolling training window and applies the same training transformation to the forecast-origin value.",
  "When X is missing, the code uses the median of the current rolling training window only; it never uses future observations.",
  "With ALLOW_X_IMPUTATION=TRUE, the external variance-regressor matrix contains standardized X and an X-missingness indicator.",
  "The Python regime-adjusted version rescales these forecasts and should be called regime-calibrated TGARCH-X, not a regime-switching TGARCH-X."
)
writeLines(
  parameterization_note,
  file.path(OUTPUT_DIR, "30_TGARCH_Parameterization_and_Naming.txt")
)

target_file <- file.path(OUTPUT_DIR, "20_Common_Target_for_Python.csv")
target_hash <- if (file.exists(target_file)) {
  digest::digest(file = target_file, algo = "sha256")
} else {
  NA_character_
}
writeLines(
  paste0("SHA-256=", target_hash),
  file.path(OUTPUT_DIR, "31_Common_Target_SHA256.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "32_R_Session_Info.txt")
)

output_files_before_manifest <- list.files(
  OUTPUT_DIR,
  full.names = TRUE,
  recursive = FALSE
)
output_manifest <- data.frame(
  File = basename(output_files_before_manifest),
  Bytes = file.info(output_files_before_manifest)$size,
  SHA256 = vapply(
    output_files_before_manifest,
    function(file_name) digest::digest(file = file_name, algo = "sha256"),
    character(1)
  ),
  stringsAsFactors = FALSE
) %>%
  arrange(File)

write_csv_safe(output_manifest, "33_Output_Manifest.csv")

# Preserve the exact sensitivity script when R exposes the --file argument.
command_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- command_arguments[grepl("^--file=", command_arguments)]
source_copy_status <- "Source script path unavailable in this execution mode."

if (length(file_argument) >= 1L) {
  source_path <- sub("^--file=", "", file_argument[[1L]])
  if (file.exists(source_path)) {
    file.copy(
      source_path,
      file.path(
        OUTPUT_DIR,
        "TGARCH_X_Reviewer_Aligned_Sensitivity_v1_00_2016_2026.R"
      ),
      overwrite = TRUE
    )
    source_copy_status <- paste0(
      "Exact source copied from: ",
      normalizePath(source_path, winslash = "/", mustWork = TRUE)
    )
  }
}

writeLines(
  c(
    "Primary R stage status: FROZEN; not replaced by this sensitivity run.",
    "Sensitivity families: WINDOW, DISTRIBUTION, X.",
    paste0("FinalTest sensitivity start: ", SENSITIVITY_FINAL_TEST_START),
    paste0("Source-copy status: ", source_copy_status)
  ),
  file.path(OUTPUT_DIR, "41_Sensitivity_Reproducibility_Note.txt")
)

# Regenerate the manifest after all sensitivity summaries and provenance files.
final_output_files <- list.files(
  OUTPUT_DIR,
  full.names = TRUE,
  recursive = FALSE
)
final_output_files <- final_output_files[
  basename(final_output_files) != "42_Final_Sensitivity_Manifest.csv"
]

final_sensitivity_manifest <- data.frame(
  File = basename(final_output_files),
  Bytes = file.info(final_output_files)$size,
  SHA256 = vapply(
    final_output_files,
    function(file_name) {
      digest::digest(file = file_name, algo = "sha256")
    },
    character(1)
  ),
  stringsAsFactors = FALSE
) %>%
  arrange(File)

write_csv_safe(
  final_sensitivity_manifest,
  "42_Final_Sensitivity_Manifest.csv"
)

# Create one clean ZIP containing all sensitivity outputs.
output_parent_directory <- normalizePath(
  dirname(OUTPUT_DIR),
  winslash = "/",
  mustWork = TRUE
)

sensitivity_zip <- file.path(
  output_parent_directory,
  paste0(
    "TGARCH_X_Reviewer_Aligned_",
    PERIOD_TAG,
    "_SENSITIVITY.zip"
  )
)

old_working_directory <- getwd()
setwd(OUTPUT_DIR)
on.exit(setwd(old_working_directory), add = TRUE)

zip_files <- list.files(
  ".",
  recursive = TRUE,
  all.files = FALSE,
  no.. = TRUE
)

if (file.exists(sensitivity_zip)) {
  file.remove(sensitivity_zip)
}

zip::zipr(
  zipfile = sensitivity_zip,
  files = zip_files,
  recurse = TRUE
)

setwd(old_working_directory)

if (!file.exists(sensitivity_zip)) {
  stop("Sensitivity ZIP was not created successfully.")
}

cat(
  "Sensitivity ZIP: ",
  normalizePath(sensitivity_zip, winslash = "/", mustWork = TRUE),
  "\n",
  sep = ""
)

elapsed_minutes <- as.numeric(
  difftime(end_time, start_time, units = "mins")
)

cat("\n=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("Run profile: ", RUN_PROFILE, "\n", sep = "")
cat("Saved forecast rows: ", nrow(results_df), "\n", sep = "")
cat("Valid GARCH forecasts: ", valid_garch_count, "\n", sep = "")
cat("Failed GARCH forecasts: ", failed_garch_count, "\n", sep = "")
cat("Runtime: ", round(elapsed_minutes, 2), " minutes\n", sep = "")
cat("QLIKE matches the Python ratio form: a/f - log(a/f) - 1.\n")
cat("Failed GARCH fits remain missing in sensitivity accuracy results.\n")
cat("Operational fallbacks are exported separately and explicitly flagged.\n")
if (exists("x_imputation_audit") && nrow(x_imputation_audit) > 0L) {
  cat(
    "TGARCH-X windows with any X imputation: ",
    sum(x_imputation_audit$N_Windows_With_Any_Imputation, na.rm = TRUE),
    " / ",
    sum(x_imputation_audit$N_Origins, na.rm = TRUE),
    "\n",
    sep = ""
  )
}
cat("Primary Python exports are intentionally not regenerated in SENSITIVITY mode.\n")
cat("Output directory: ", OUTPUT_DIR, "\n", sep = "")
cat("=============================================================================\n")
