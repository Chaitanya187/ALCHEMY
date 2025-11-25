"""
Phase 1: Data Acquisition & Preprocessing
Fetch compound data from ChEMBL and preprocess SMILES strings
"""

import os
import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import time

class ChEMBLDataFetcher:
    """Fetch active compounds for a specific target from ChEMBL"""
    
    def __init__(self, target_name, activity_threshold=10000):
        """
        Initialize the fetcher
        
        Args:
            target_name: Name of the target protein (e.g., 'EGFR', 'CDK2')
            activity_threshold: IC50/Ki threshold in nM (default: 10000 nM = 10 µM)
        """
        self.target_name = target_name
        self.activity_threshold = activity_threshold
        self.target_client = new_client.target
        self.activity_client = new_client.activity
        
    def find_target_chembl_id(self):
        """Search for target and return ChEMBL ID"""
        print(f"Searching for target: {self.target_name}...")
        targets = self.target_client.search(self.target_name)
        
        if not targets:
            print(f"No targets found for {self.target_name}")
            return None
            
        # Get the first human target
        for target in targets:
            if target.get('organism', '').lower() == 'homo sapiens':
                chembl_id = target['target_chembl_id']
                pref_name = target.get('pref_name', 'Unknown')
                print(f"Found target: {pref_name} ({chembl_id})")
                return chembl_id
        
        # If no human target, take the first one
        chembl_id = targets[0]['target_chembl_id']
        pref_name = targets[0].get('pref_name', 'Unknown')
        print(f"Found target: {pref_name} ({chembl_id})")
        return chembl_id
    
    def fetch_bioactivity_data(self, target_chembl_id, limit=10000):
        """Fetch bioactivity data for the target"""
        print(f"Fetching bioactivity data for {target_chembl_id}...")
        
        try:
            # Try to get IC50 data with limit
            activities = self.activity_client.filter(
                target_chembl_id=target_chembl_id,
                standard_type='IC50'
            ).only(['molecule_chembl_id', 'canonical_smiles', 'standard_value', 
                    'standard_units', 'standard_type'])
            
            print("Processing activities...")
            data = []
            count = 0
            
            for act in activities:
                if count >= limit:
                    break
                    
                if act.get('canonical_smiles') and act.get('standard_value'):
                    try:
                        value = float(act['standard_value'])
                        units = act.get('standard_units', 'nM')
                        
                        # Convert to nM if needed
                        if units == 'uM':
                            value = value * 1000
                        elif units == 'pM':
                            value = value / 1000
                        
                        # Filter by activity threshold
                        if value <= self.activity_threshold:
                            data.append({
                                'molecule_chembl_id': act['molecule_chembl_id'],
                                'canonical_smiles': act['canonical_smiles'],
                                'IC50_nM': value
                            })
                    except (ValueError, TypeError):
                        continue
                
                count += 1
                if count % 100 == 0:
                    print(f"  Processed {count} activities, found {len(data)} active compounds...")
            
            df = pd.DataFrame(data)
            print(f"Found {len(df)} active compounds (IC50 <= {self.activity_threshold/1000} µM)")
            return df
            
        except Exception as e:
            print(f"Error fetching data: {e}")
            return pd.DataFrame()
    
    def fetch_data(self):
        """Main method to fetch data"""
        target_id = self.find_target_chembl_id()
        if not target_id:
            return None
        
        df = self.fetch_bioactivity_data(target_id)
        return df


class SMILESPreprocessor:
    """Clean and filter SMILES strings"""
    
    def __init__(self, min_mw=180, max_mw=500):
        """
        Initialize preprocessor
        
        Args:
            min_mw: Minimum molecular weight (default: 180)
            max_mw: Maximum molecular weight (default: 500)
        """
        self.min_mw = min_mw
        self.max_mw = max_mw
        self.allowed_elements = {'C', 'H', 'O', 'N', 'S', 'F', 'Cl', 'Br'}
        self.uncharger = rdMolStandardize.Uncharger()
        
    def canonicalize_smiles(self, smiles):
        """Convert SMILES to canonical form"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol, canonical=True)
        except:
            return None
    
    def remove_salts_and_standardize(self, smiles):
        """Remove salts and standardize molecule"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            # Remove salts (keep largest fragment)
            remover = rdMolStandardize.LargestFragmentChooser()
            mol = remover.choose(mol)
            
            # Remove charges
            mol = self.uncharger.uncharge(mol)
            
            # Remove stereochemistry
            Chem.RemoveStereochemistry(mol)
            
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        except Exception as e:
            return None
    
    def check_element_composition(self, smiles):
        """Check if molecule contains only allowed elements"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            
            atoms = set([atom.GetSymbol() for atom in mol.GetAtoms()])
            return atoms.issubset(self.allowed_elements)
        except:
            return False
    
    def passes_lipinski(self, smiles):
        """Check Lipinski's Rule of 5"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            # Lipinski's Rule of 5
            if mw < self.min_mw or mw > self.max_mw:
                return False
            if logp > 5:
                return False
            if hbd > 5:
                return False
            if hba > 10:
                return False
            
            return True
        except:
            return False
    
    def is_valid_molecule(self, smiles):
        """Check if SMILES is valid and passes all filters"""
        if not smiles or pd.isna(smiles):
            return False
        
        # Try to parse
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # Check element composition
        if not self.check_element_composition(smiles):
            return False
        
        # Check Lipinski's Rule
        if not self.passes_lipinski(smiles):
            return False
        
        return True
    
    def process_dataframe(self, df):
        """Process entire dataframe of SMILES"""
        print("\n=== Starting SMILES Preprocessing ===")
        print(f"Initial compounds: {len(df)}")
        
        # Remove duplicates
        df = df.drop_duplicates(subset=['canonical_smiles'])
        print(f"After removing duplicates: {len(df)}")
        
        # Standardize and clean
        print("Standardizing molecules (removing salts, stereochemistry)...")
        df['clean_smiles'] = df['canonical_smiles'].apply(self.remove_salts_and_standardize)
        df = df.dropna(subset=['clean_smiles'])
        print(f"After standardization: {len(df)}")
        
        # Remove duplicates again (some might become identical after standardization)
        df = df.drop_duplicates(subset=['clean_smiles'])
        print(f"After removing duplicate standardized structures: {len(df)}")
        
        # Filter by validity
        print("Applying filters (elements, Lipinski's Rule, MW)...")
        df['is_valid'] = df['clean_smiles'].apply(self.is_valid_molecule)
        df = df[df['is_valid'] == True]
        print(f"After applying all filters: {len(df)}")
        
        return df


def main():
    """Main execution function"""
    print("="*60)
    print("PHASE 1: DATA ACQUISITION & PREPROCESSING")
    print("="*60)
    
    # List of targets to try (in order of preference)
    targets = [
        'Epidermal growth factor receptor erbB1',
        'Cyclin-dependent kinase 2', 
        'Dopamine D2 receptor',
        'Acetylcholinesterase',
        'Estrogen receptor alpha',
        'Cannabinoid CB1 receptor'
    ]
    min_required_compounds = 2000
    
    df_final = None
    target_used = None
    
    for target in targets:
        print(f"\n{'='*60}")
        print(f"Attempting to fetch data for: {target}")
        print(f"{'='*60}")
        
        # Fetch data
        fetcher = ChEMBLDataFetcher(target, activity_threshold=10000)
        df = fetcher.fetch_data()
        
        if df is None or len(df) == 0:
            print(f"❌ No data found for {target}. Trying next target...")
            continue
        
        # Preprocess
        preprocessor = SMILESPreprocessor(min_mw=180, max_mw=500)
        df_clean = preprocessor.process_dataframe(df)
        
        if len(df_clean) >= min_required_compounds:
            print(f"✓ SUCCESS: Found {len(df_clean)} valid compounds for {target}")
            df_final = df_clean
            target_used = target
            break
        else:
            print(f"⚠ Warning: Only {len(df_clean)} compounds found for {target} (need {min_required_compounds})")
            if df_final is None or len(df_clean) > len(df_final):
                df_final = df_clean
                target_used = target
    
    if df_final is None or len(df_final) == 0:
        print("\n❌ ERROR: Could not fetch sufficient data from any target!")
        return
    
    # Save results
    print(f"\n{'='*60}")
    print("SAVING RESULTS")
    print(f"{'='*60}")
    
    # Save clean SMILES to text file
    output_file = '../data/clean_smiles.txt'
    with open(output_file, 'w') as f:
        for smiles in df_final['clean_smiles']:
            f.write(smiles + '\n')
    
    print(f"✓ Saved {len(df_final)} SMILES to '{output_file}'")
    
    # Save detailed CSV
    csv_file = '../data/processed_compounds.csv'
    df_final[['molecule_chembl_id', 'canonical_smiles', 'clean_smiles', 'IC50_nM']].to_csv(
        csv_file, index=False
    )
    print(f"✓ Saved detailed data to '{csv_file}'")
    
    # Calculate and display statistics
    print(f"\n{'='*60}")
    print("DATASET STATISTICS")
    print(f"{'='*60}")
    print(f"Target Used: {target_used}")
    print(f"Total Clean SMILES: {len(df_final)}")
    print(f"Unique Compounds: {df_final['clean_smiles'].nunique()}")
    
    # Calculate molecular properties for statistics
    print("\nCalculating molecular properties...")
    mws = []
    logps = []
    
    for smiles in df_final['clean_smiles'].head(100):  # Sample for speed
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mws.append(Descriptors.MolWt(mol))
                logps.append(Descriptors.MolLogP(mol))
        except:
            continue
    
    if mws:
        print(f"\nMolecular Weight: {min(mws):.1f} - {max(mws):.1f} (avg: {sum(mws)/len(mws):.1f})")
    if logps:
        print(f"LogP: {min(logps):.2f} - {max(logps):.2f} (avg: {sum(logps)/len(logps):.2f})")
    
    print(f"\n{'='*60}")
    print("✓ PHASE 1 COMPLETE!")
    print(f"{'='*60}")
    print(f"\nOutput files:")
    print(f"  - {output_file} (for LSTM training)")
    print(f"  - {csv_file} (detailed compound data)")


if __name__ == "__main__":
    main()
