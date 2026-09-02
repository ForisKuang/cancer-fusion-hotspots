from datetime import datetime
from typing import Optional

from pydantic import BaseModel, ConfigDict


class AlgorithmResult(BaseModel):
    """Structured output of a single hotspot-detection algorithm run."""

    model_config = ConfigDict(extra="forbid")

    Algorithm: str
    Algorithm_version: Optional[str] = None
    Parameters: Optional[dict] = None
    Summary: Optional[dict] = None
    Tables: Optional[dict] = None
    Warnings: Optional[list] = None
    Created_at: Optional[datetime] = None
    Input_fingerprint: Optional[str] = None
