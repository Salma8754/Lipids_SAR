from sqlalchemy import Column, String, Float

from data_client import CustomBase


class EpoMerged(CustomBase):
    __tablename__ = "VDP7_epo_merged"
    __table_args__ = {"schema": "PSA_RESULTS_MRNA"}

    study_id = Column(String)
    study_group_id = Column(String)
    catlipid_id = Column(String, primary_key=True)
    peg_lipid = Column(String)
    cholesterol_lipid = Column(String)
    helper_lipid = Column(String)
    composition_cat = Column(Float)
    composition_help = Column(Float)
    composition_chol = Column(Float)
    composition_peg = Column(Float)
    roa = Column(String)
    species = Column(String)
    data_owner = Column(String)
    epo_conc_mean_6hr = Column(Float)
    epo_conc_mean_24hr = Column(Float)
    experiment_count = Column(Float)
    pdi = Column(Float)
    encapsulation = Column(Float)
    size = Column(Float)
    mixing = Column(String)
    formulation_type = Column(String)
    final_buffer = Column(String)
    buffer_exchange = Column(String)
    formulation_buffer = Column(String)
    dose = Column(Float)
    n_p_ratio = Column(Float)
    mrna = Column(String)
    conc = Column(Float)
    structure_smiles = Column(String)
    cat_lipid_nm = Column(String)
    canonical_smiles = Column(String)


class catlipids(CustomBase):
    __tablename__ = "VDP7_catlipids_merged"

    catlipid_id = Column(String, primary_key=True)
    mol_wt = Column(String)
    catlipid_family_name = Column(String)

