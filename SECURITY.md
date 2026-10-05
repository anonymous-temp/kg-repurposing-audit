# Security and data handling

Never report a security issue by posting credentials or restricted records. Provide a minimal synthetic example and the affected release. For a private disclosure, contact the corresponding author at wangzeyuan@sentumhealth.com.

This repository has no hosted patient-data service. Its default workflow reads caller-owned local inputs and writes caller-selected output directories. Database folders, checkpoints and outputs are excluded from Git. Users remain responsible for their local data permissions and storage controls.

Do not load untrusted model checkpoints. The pipeline writes PyTorch state dictionaries, not a public collection of clinical models. No secrets should be required for the synthetic tests or public registry API.
