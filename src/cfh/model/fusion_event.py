from typing import Optional

from pydantic import BaseModel, ConfigDict


class FusionEvent(BaseModel):
    """A single normalized structural-variant / fusion event for one sample."""

    model_config = ConfigDict(extra="forbid")

    Event_id: str
    Cohort: str
    Sequencing_panel_id: Optional[str] = None
    Sample_id: Optional[str] = None
    Patient_id: Optional[str] = None
    Site1_gene: Optional[str] = None
    Site2_gene: Optional[str] = None
    Five_prime_gene: Optional[str] = None
    Three_prime_gene: Optional[str] = None
    Fusion_name: Optional[str] = None
    Event_class: Optional[str] = None
    Connection_type: Optional[str] = None
    Frame_status: Optional[str] = None
    Is_protein_fusion: Optional[bool] = None
    Is_antisense: Optional[bool] = None
    Confidence_class: Optional[str] = None
    Paired_end_read_support: Optional[int] = None
    Split_read_support: Optional[int] = None
    Tumor_variant_count: Optional[int] = None
    Site1_description: Optional[str] = None
    Site2_description: Optional[str] = None
    Annotation: Optional[str] = None
    Event_info: Optional[str] = None
    Source_row_number: Optional[int] = None
