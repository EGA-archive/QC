from multiqc.base_module import BaseMultiqcModule
import logging
import json
import os
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Median Coverage Across Reference",
            anchor="genome-reference",
            info="This plot displays sequencing coverage across chromosomes, summarized as the median coverage within non-overlapping 3 Mb windows. It helps identify regions with unusually high or low coverage, which may reflect technical biases, duplications, or low-complexity regions"
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

        coverage_median_data = {}

        # === Parse coverage_median y create section ===
        log.info("coverage_median_3Mb_complete.json...")
        # === Parsear archivo de cabecera BAM ===
        for f in self.find_log_files("genomewide_coverage"):
            log.info(f"Archivo encontrado: {f['fn']}")
            coverage_median_data = json.dumps(f["f"])

        log.info(f"Parsed coverage_median data")    

        self.add_section( 
            anchor="median_coverage_across_reference",
            content=html_import(
                 os.path.join(
                    os.path.dirname(__file__), "assets", "html", "coverage.html"
                ),
                """JSON.parse("{{ coverage_intervals }}")""",
                f"""JSON.parse({coverage_median_data})"""
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
