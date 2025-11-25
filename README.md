# De Novo Drug Discovery using LSTM

**Target:** EGFR (Epidermal Growth Factor Receptor)  
**Approach:** Character-level LSTM for SMILES generation  
**Dataset:** 2,664 drug-like molecules from ChEMBL

---

## 📁 Project Structure

```
DENOVO/
│
├── README.md                  # This file - Project overview
├── Plan.md                    # Complete project plan & timeline
├── requirements.txt           # Python dependencies
├── .gitignore                 # Git ignore rules
│
├── data/                      # Training & output data
│   ├── README.md              # Data directory documentation
│   ├── clean_smiles.txt       # 2,664 SMILES for LSTM training
│   └── processed_compounds.csv # Metadata & ChEMBL references
│
├── scripts/                   # Python scripts
│   ├── README.md              # Scripts documentation
│   ├── phase1_data_acquisition.py  # ChEMBL data fetching & preprocessing
│   ├── fetch.py               # Original fetch script
│   ├── phase2_train_lstm.py   # LSTM model training (TODO)
│   ├── phase3_generate.py     # Molecule generation (TODO)
│   └── phase4_evaluate.py     # Filtering & evaluation (TODO)
│
├── models/                    # Saved model weights
│   └── README.md              # Models directory documentation
│
├── results/                   # Generated molecules & analysis
│   ├── README.md              # Results documentation
│   ├── generated_molecules/   # Generated SMILES files
│   └── plots/                 # Visualizations
│
└── docs/                      # Documentation
    ├── README.md              # Documentation index
    ├── PHASE1_SUMMARY.md      # Detailed Phase 1 report
    ├── PHASE1_FLOW.md         # Visual flow diagrams
    ├── PHASE1_QUICK_REFERENCE.md  # Quick reference guide
    └── WORKFLOW.md            # Additional workflow notes
```

---

## 🎯 Project Status

### Phase 1: Data Acquisition & Preprocessing ✅ COMPLETE
- **Target:** EGFR (CHEMBL203)
- **Dataset:** 2,664 drug-like molecules
- **Quality:** 100% valid, Lipinski-compliant
- **Files:** `data/clean_smiles.txt`, `data/processed_compounds.csv`

### Phase 2: LSTM Model Training 🔄 TODO
- Character vocabulary extraction
- LSTM architecture (2 layers, 256-512 units)
- Training for 50-100 epochs

### Phase 3: Molecule Generation 🔄 TODO
- Generate 10,000 novel SMILES
- Temperature sampling (0.7-1.0)

### Phase 4: Evaluation & Filtering 🔄 TODO
- Validity check (RDKit)
- Novelty assessment
- Drug-likeness scoring (QED, SA Score)

---

## 🚀 Quick Start

### Phase 1 (Completed)
```bash
cd C:\DENOVO
python scripts/phase1_data_acquisition.py
```

### Phase 2 (Next Step)
```bash
# Train LSTM model
python scripts/phase2_train_lstm.py
```

---

## 📊 Dataset Statistics

| Metric | Value |
|--------|-------|
| Target Protein | EGFR (Epidermal Growth Factor Receptor) |
| Total Molecules | 2,664 |
| Molecular Weight | 202-498 Da (avg: 360 Da) |
| LogP | 0.49-4.77 (avg: 3.43) |
| Activity (IC50) | All ≤ 10 µM |
| Drug-likeness | 100% Lipinski-compliant |

---

## 📚 Documentation

### Main Documentation
- **[Plan.md](Plan.md)** - Complete project roadmap with weekly schedule
- **[requirements.txt](requirements.txt)** - Python package dependencies

### Phase 1 Documentation (Completed)
- **[docs/PHASE1_SUMMARY.md](docs/PHASE1_SUMMARY.md)** - Comprehensive Phase 1 report
- **[docs/PHASE1_FLOW.md](docs/PHASE1_FLOW.md)** - Visual workflow diagrams
- **[docs/PHASE1_QUICK_REFERENCE.md](docs/PHASE1_QUICK_REFERENCE.md)** - Quick reference guide
- **[docs/WORKFLOW.md](docs/WORKFLOW.md)** - Additional workflow notes

### Directory Documentation
- **[data/README.md](data/README.md)** - Training data and files explanation
- **[scripts/README.md](scripts/README.md)** - Scripts overview and usage
- **[models/README.md](models/README.md)** - Model storage structure
- **[results/README.md](results/README.md)** - Results and outputs explanation
- **[docs/README.md](docs/README.md)** - Documentation index

---

## 🛠️ Requirements

```bash
# Install all dependencies
pip install -r requirements.txt

# Or install individually:
# Core packages (Phase 1)
pip install chembl_webresource_client rdkit pandas numpy

# Deep Learning (Phase 2+)
pip install torch matplotlib seaborn
```

---

## 📖 References

- ChEMBL Database: https://www.ebi.ac.uk/chembl/
- RDKit Documentation: https://www.rdkit.org/docs/
- SMILES Notation: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system

---

## 👤 Author

**Student Project:** De Novo Drug Discovery  
**Institution:** [Your University]  
**Timeline:** 8-9 Weeks (Semester Project)
