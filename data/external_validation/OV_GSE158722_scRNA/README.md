# OV_GSE158722 scRNA source matrix

`OV_GSE158722_expression.h5` is stored as split parts because the original
file is larger than GitHub's single Git LFS object limit.

To restore the original source file on Windows PowerShell:

```powershell
Get-Content -Encoding Byte -Path `
  .\OV_GSE158722_expression.h5.part001,`
  .\OV_GSE158722_expression.h5.part002 |
  Set-Content -Encoding Byte .\OV_GSE158722_expression.h5
```

Expected SHA256 for the restored file:

```text
d9dfe861d9dff40abb96bc971e163e1500487039c4921f186bab77651c503860
```

