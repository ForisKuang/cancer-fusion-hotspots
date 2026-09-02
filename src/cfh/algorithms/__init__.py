"""Built-in hotspot-detection algorithm plugins."""

# Import built-ins so their registration decorators run on package import.
from cfh.algorithms.joint_partner import JointPartnerMode

__all__ = ["JointPartnerMode"]
