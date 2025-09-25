from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table
import logging
import re

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Summary Tables",
            anchor="summary-tables",
        )

        header_data  = {}
        genome_data  = {}
        sample_name  = None

        # 1) Parseamos bam_header_info.txt (TableMaker / parse_tabbed_key_value)
        for f in self.find_log_files("tablemaker/header"):
            log.info("Parsing BAM header file")
            header_data = parse_tabbed_key_value(f["f"])
            # Extraemos el valor de 'Sample' (p.ej. "HG00096") y lo almacenamos en sample_name
            # Al mismo tiempo, lo eliminamos del diccionario para que no aparezca como columna
            sample_name = header_data.pop("Sample", None)
            if not sample_name:
                sample_name = "Sample"  # fallback, por si no hubiera campo Sample

        if sample_name is None:
            # Por si no encontramos nigún bam_header_info.txt con campo "Sample"
            sample_name = "Sample"

        log.info("Buscando genome_results.txt...")
        for f in self.find_log_files("tablemaker/genome"):
            log.info(f"Archivo encontrado: {f['fn']}")
            genome_data = parse_tabbed_key_value(f["f"])

        log.info(f"Parsed genome data: {genome_data}")
        log.info(f"Usando sample_name = {sample_name} para ambas tablas")

        # 2) Construimos la primera tabla ("BAM Header Info"), usando sample_name como índice
        header_table_data = { sample_name: header_data }
        # Generamos los metadatos de cabeceras de columnas (Título + rid) con las claves de header_data
        header_headers = {
            k: { "title": k, "rid": f"bam_header_table-{k.replace(' ', '_')}" }
            for k in header_data.keys()
        }

        self.add_section(
            name="BAM Header Info",
            anchor="tablemaker_bam",
            description="Basic metadata extracted from the BAM header.",
            plot=table.plot(
                header_table_data,
                headers=header_headers,
                pconfig={
                    "namespace": "bam_header_table",
                    "id": "bam_header",
                    "title": "bam_header"
                }
            )
        )

        # 3) Construimos la segunda tabla ("Genome Results Summary"), reutilizando el mismo sample_name
        genome_table_data = { sample_name: genome_data }
        genome_headers = {
            k: { "title": k, "rid": f"genome_results_table-{k.replace(' ', '_')}" }
            for k in genome_data.keys()
        }

        self.add_section(
            name="Genome Results Summary",
            anchor="tablemaker_genome",
            description="General information from the BAM file.",
            plot=table.plot(
                genome_table_data,
                headers=genome_headers,
                pconfig={
                    "namespace": "genome_results_table",
                    "id": "genome_results",
                    "title": "genome_results"
                }
            )
        )

def parse_tabbed_key_value(file_contents):
    data = {}
    for line in file_contents.strip().splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 2:
            key, value = parts
            data[key.strip()] = value.strip()
    return data
