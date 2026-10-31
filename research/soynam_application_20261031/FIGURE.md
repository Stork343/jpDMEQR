# Empirical result figure

The figure asks how canopy coverage relates to yield after genomic adjustment and whether it improves prediction for held-out lines. Panel A displays the six year-specific coefficients. Panel B displays the three registered comparisons of the joint predictor with and without canopy. Its simultaneous intervals retain the error allocation across all six planned comparisons.

| Element | Meaning | Source |
| --- | --- | --- |
| Panel A points | Complete profile coefficients in bu/ac per ten percentage points of canopy coverage | `results/tables/soynam_20261031/inference.csv`, reference setting |
| Panel A intervals | Bonferroni simultaneous intervals from complete-line influences, family size eleven | Saved joint coordinate and contrast covariance |
| Panel B points | Percentage reduction in calibrated held-out pinball loss when canopy is added | `results/tables/soynam_20261031/prediction_comparisons.csv`, calibrated joint-genomic comparison |
| Panel B intervals | Paired-line delta-method intervals, family size six | Same comparison table |

The inference sample has 5533 independent line clusters under the conditional field model. Prediction uses the same 1107 test lines for every method and conditions on the fitted training models. No row or result is omitted from either displayed family. Plotting uses already aggregated values and does not estimate additional intervals.

Colour and marker shape both identify the year or comparator. The coefficient axis includes zero, and the prediction panel displays the zero-improvement reference. The two panels have different units. Their axes explicitly state those units. `plot_results.R` produces the editable R plot objects, a vector PDF for the manuscript and a 300 dpi PNG.

Verification compares plotted row counts with the complete table families. PDF checks use extracted text and page coordinates. Visual inspection is excluded by the user's instruction.
