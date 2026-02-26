"""
DockIt — Email notification module.

Sends email notifications when long-running jobs (e.g. deep screening) complete.
Requires SMTP environment variables to be configured. If not configured,
notifications are silently skipped.
"""

from __future__ import annotations

import logging
import os
import smtplib
from email.mime.text import MIMEText

logger = logging.getLogger(__name__)

SMTP_HOST: str = os.environ.get("SMTP_HOST", "")
SMTP_PORT: int = int(os.environ.get("SMTP_PORT", "587"))
SMTP_USER: str = os.environ.get("SMTP_USER", "")
SMTP_PASS: str = os.environ.get("SMTP_PASS", "")
FROM_EMAIL: str = os.environ.get("FROM_EMAIL", "noreply@dockit.app")


def send_notification_email(to_email: str, job_id: str, results_url: str = "") -> bool:
    """Send email notification when a screening job completes.

    Parameters
    ----------
    to_email : str
        Recipient email address.
    job_id : str
        The job UUID.
    results_url : str, optional
        Direct URL to the results page. Defaults to localhost URL.

    Returns
    -------
    bool
        True if the email was sent successfully, False otherwise.
    """
    if not SMTP_HOST or not to_email:
        logger.info("Email notification skipped (no SMTP config or no email)")
        return False

    subject = f"DockIt - Your screening is complete (Job {job_id[:8]})"
    body = f"""Hello,

Your DockIt molecular screening is complete.

Job ID: {job_id}
Results: {results_url or f'http://localhost:3000?job={job_id}'}

View the results in the DockIt interface.

Best regards,
The DockIt Team
"""
    try:
        msg = MIMEText(body, "plain", "utf-8")
        msg["Subject"] = subject
        msg["From"] = FROM_EMAIL
        msg["To"] = to_email

        with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
            if SMTP_PORT == 587:
                server.starttls()
            if SMTP_USER:
                server.login(SMTP_USER, SMTP_PASS)
            server.sendmail(FROM_EMAIL, [to_email], msg.as_string())

        logger.info("Notification email sent to %s for job %s", to_email, job_id)
        return True
    except Exception as exc:
        logger.warning("Failed to send notification email: %s", exc)
        return False
