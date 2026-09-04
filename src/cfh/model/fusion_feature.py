from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


class DomainRetentionDetail(BaseModel):
    """Quantitative retention of one inclusive protein-domain interval."""

    model_config = ConfigDict(extra="forbid")

    Domain_start_aa: Optional[int] = None
    Domain_end_aa: Optional[int] = None
    Retained_start_aa: Optional[int] = None
    Retained_end_aa: Optional[int] = None
    Retained_fraction: Optional[float] = Field(default=None, ge=0.0, le=1.0)
    Is_truncated: Optional[bool] = None


class FusionFeature(BaseModel):
    """Transcript/exon/domain-level annotation for one gene's role in a FusionEvent."""

    model_config = ConfigDict(extra="forbid")

    Event_id: str
    Gene: str
    Role: Optional[str] = None
    Transcript_id: Optional[str] = None
    Breakpoint_exon: Optional[int] = None
    Breakpoint_intron: Optional[int] = None
    Retained_exons: Optional[list] = None
    Lost_exons: Optional[list] = None
    Junction_position_aa: Optional[int] = None
    Retained_domains: Optional[list] = None
    Lost_domains: Optional[list] = None
    Disrupted_domains: Optional[list] = None
    Domain_retention_flags: Optional[dict] = None
    Domain_retention_details: Optional[dict[str, DomainRetentionDetail]] = None
