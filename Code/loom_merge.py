# merge CKG and CKN into one loom file
import loompy
files=['/mnt/c/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/CKG.loom','/mnt/c/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/CKN.loom']
output_filename='/mnt/c/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/scVelo/CK.loom'
loompy.combine(files, output_filename, key="Accession")
