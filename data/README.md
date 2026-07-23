# Bundled IHS classification

`ihs_classification.bin` is a compact, lossless representation of the labeled
616 by 616 `EssMatIHS.csv` supplied with this project. It stores 56 unique
customer-industry profiles plus their column map; `ihs_industry_codes.txt`
stores the original row/column order. The source CSV is byte-identical to
`ch-diem/esri_tutorial/data/ess_mat_n4_ihs.csv` at SHA-256
`10db181089d8fe3a528be7e03c0b1ab2b13dbe72da6cf627caeb16f76ceae014`.

The package exposes these assets through `ihs_input_classification()` and
`ihs_industry_codes()`. The existing `IndustryInfo(industry_of_firm)` default
remains all-essential for backward compatibility.
