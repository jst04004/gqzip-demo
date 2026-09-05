"""
GQZip Cryptographic License Verification and Local Metering Engine.
"""

import os
import sys
import json
import base64
import hmac
import hashlib
import time
from typing import Dict, Tuple, Optional

# Secret HMAC seed for signature verification
_LICENSE_SECRET = b"GQZip-Master-Secret-Key-2026-Production"
DEFAULT_ALLOWANCE_BYTES = 1024 * 1024 * 1024 * 1024  # 1 Terabyte (1 TB)

def get_metering_dir() -> str:
    """Returns persistent platform-appropriate metering directory."""
    if os.name == "nt":
        base = os.environ.get("USERPROFILE", os.path.expanduser("~"))
    else:
        base = os.path.expanduser("~")
    path = os.path.join(base, ".gqzip")
    os.makedirs(path, exist_ok=True)
    return path

def get_metering_filepath() -> str:
    return os.path.join(get_metering_dir(), "metering.json")

def get_license_filepath() -> str:
    return os.path.join(get_metering_dir(), "license.key")

def load_metering_state() -> Dict:
    """Loads persistent local byte metering state."""
    fp = get_metering_filepath()
    if os.path.exists(fp):
        try:
            with open(fp, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass
    return {
        "total_bytes_compressed": 0,
        "allowance_limit_bytes": DEFAULT_ALLOWANCE_BYTES,
        "installation_timestamp": time.time()
    }

def save_metering_state(state: Dict) -> None:
    """Saves persistent local byte metering state."""
    fp = get_metering_filepath()
    try:
        with open(fp, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=2)
    except Exception:
        pass

def update_processed_bytes(added_bytes: int) -> Dict:
    """Increments total compressed bytes counter."""
    state = load_metering_state()
    state["total_bytes_compressed"] = state.get("total_bytes_compressed", 0) + added_bytes
    save_metering_state(state)
    return state

def verify_license_key(key_str: str) -> Tuple[bool, str, Dict]:
    """
    Verifies cryptographic signature, payload structure, and expiration of a license key.
    Returns: (is_valid, message, payload_dict)
    """
    key_str = key_str.strip()
    if not key_str.startswith("GQZIP-v1-"):
        return False, "Invalid license key format.", {}
    
    raw = key_str[len("GQZIP-v1-"):]
    if "." not in raw:
        return False, "Invalid license signature structure.", {}
    
    payload_b64, signature = raw.split(".", 1)
    
    # Verify HMAC-SHA256 signature
    expected_sig = hmac.new(_LICENSE_SECRET, payload_b64.encode("ascii"), hashlib.sha256).hexdigest()
    if not hmac.compare_digest(expected_sig, signature):
        return False, "Cryptographic signature verification failed.", {}
    
    # Decode JSON payload
    try:
        payload_bytes = base64.urlsafe_b64decode(payload_b64.encode("ascii"))
        payload = json.loads(payload_bytes.decode("utf-8"))
    except Exception as e:
        return False, f"Failed to parse license payload: {e}", {}
    
    # Check expiration date if specified
    exp_str = payload.get("expiry")
    if exp_str:
        try:
            exp_time = time.mktime(time.strptime(exp_str, "%Y-%m-%d"))
            if time.time() > exp_time:
                return False, f"License key expired on {exp_str}.", payload
        except Exception:
            pass
            
    org = payload.get("org", "Authorized User")
    l_type = payload.get("type", "standard").upper()
    return True, f"Valid {l_type} License issued to {org}.", payload

def get_active_license() -> Tuple[bool, str, Optional[Dict]]:
    """Checks environment variable or local license key file."""
    env_key = os.environ.get("GQZIP_LICENSE_KEY")
    if env_key:
        valid, msg, payload = verify_license_key(env_key)
        if valid:
            return True, f"[ENV] {msg}", payload
    
    license_fp = get_license_filepath()
    if os.path.exists(license_fp):
        try:
            with open(license_fp, "r", encoding="utf-8") as f:
                file_key = f.read().strip()
            valid, msg, payload = verify_license_key(file_key)
            if valid:
                return True, f"[FILE] {msg}", payload
        except Exception:
            pass
            
    return False, "No active license key found (Universal 1 TB Freemium Allowance Active).", None

def check_allowance_permission(incoming_bytes: int = 0) -> Tuple[bool, str]:
    """
    Checks if incoming compression job is permitted under 1 TB Freemium Allowance or Active License.
    """
    has_license, lic_msg, payload = get_active_license()
    if has_license:
        return True, f"Licensed: {lic_msg}"
        
    state = load_metering_state()
    current = state.get("total_bytes_compressed", 0)
    limit = state.get("allowance_limit_bytes", DEFAULT_ALLOWANCE_BYTES)
    
    if current + incoming_bytes > limit:
        msg = (
            f"\n========================================================================================\n"
            f" [GQZIP FREEMIUM NOTICE] You have processed {current / (1024**4):.2f} TB of raw FASTQ data.\n"
            f" Your Universal 1 Terabyte (1 TB) Free Evaluation Allowance is complete!\n\n"
            f" * ACADEMIC / NON-PROFIT LABS: Request your free annual renewal at contact@gqzip.org\n"
            f" * ENTERPRISE / COMMERCIAL LABS: Activate your commercial license key with `gqzip --license KEY`\n"
            f"========================================================================================\n"
        )
        return False, msg
        
    remaining_gb = (limit - current) / (1024**3)
    return True, f"Freemium Allowance Active ({remaining_gb:.1f} GB remaining of 1 TB allowance)"

def save_license_key(key_str: str) -> Tuple[bool, str]:
    """Saves valid license key to local user config."""
    valid, msg, payload = verify_license_key(key_str)
    if not valid:
        return False, msg
    
    license_fp = get_license_filepath()
    with open(license_fp, "w", encoding="utf-8") as f:
        f.write(key_str.strip())
        
    return True, f"Successfully activated license key for {payload.get('org', 'User')}!"
