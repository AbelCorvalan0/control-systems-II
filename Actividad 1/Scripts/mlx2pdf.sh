# Run mlx2pdf for all files.
#matlab -batch "convertMLXFiles2pdf"

# Specific file convertion
matlab -batch "matlab.internal.liveeditor.openAndConvert('actividad1_4b.mlx', 'output.pdf'); exit"


