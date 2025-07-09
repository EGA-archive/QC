from multiqc.base_module import BaseMultiqcModule
import logging
import json
import os
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Median Mapping Quality Across Reference",
            anchor="genome-reference",
            info="This plot shows mapping quality across chromosomes, summarized as the median mapping quality within non‐overlapping 3 Mb windows. Regions with consistently low mapping quality may point to repetitive or ambiguous sequence contexts, while areas of high quality indicate confidently aligned reads"
        )

        log.info("preparing static files...")
        # add static files
        self.css = {
            "assets/css/genomewide.css": os.path.join(
                os.path.dirname(__file__), "assets", "css", "genomewide.css"
            )
        }
        self.js = {
            "assets/js/egafile-curve-multiple-overlaping-scrollable.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "egafile-curve-multiple-overlaping-scrollable.js"
            ),
            "assets/js/egafile-curve-multiple-overlaping.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "egafile-curve-multiple-overlaping.js"
            )
        }

        mapq_median_data = {}
        
        # === Parse mapq_median y create section ===
        log.info("mapq_median_3Mb_complete.json...")
        for f in self.find_log_files("genomewide_mapq"):
            log.info(f"Archivo encontrado: {f['fn']}")
            mapq_median_data = json.dumps(f["f"])

        log.info(f"Parsed mapq_median data")


        self.add_section(
            anchor="median_mapping_quality_across_reference",
            content=html_import(
                os.path.join(
                    os.path.dirname(__file__), "assets", "html", "mapq.html"
                ),
                """JSON.parse("{{ mapq_intervals }}")""",
                f"""JSON.parse({mapq_median_data})"""
            )
        )


#imports html from a file and replaces all occurrences of input_substring to output_substring
def html_import(filename, input_substring , output_substring):
    log.info(filename)
    # Open the file in read mode
    file = open(filename, "r")
    
    # Read the entire content of the file
    content = file.read()
    content = content.replace(input_substring, output_substring)    
    # Close the file
    file.close()
    # Return content
    return content        
