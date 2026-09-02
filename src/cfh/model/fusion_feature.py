from typing import Optional

from pydantic import BaseModel, ConfigDict


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
