Trabajar directamente con CRAM es complicado porque muchas herramientas no lo soportan del todo bien. Por eso decidimos convertir los CRAM a BAM con samtools view -b, que permite hacer esta conversión siempre que tengas acceso al MD5sum de la referencia. Esto es obligatorio en CRAM a menos que la secuencia FASTA esté embebida en el archivo. Herramientas que usamos: 

Qualimap:
 Saca estadísticas muy completas sobre los BAM, como calidad de mapeo, cobertura, composición de nucleótidos, presencia de homopolímeros, etc. También genera un archivo resumen (genome_results.txt) con métricas generales.


RSeQC:
 Crea gráficos sobre dónde caen las lecturas en el genoma (exones, intrones, UTRs…), usando un archivo BED con las anotaciones.


Samtools (depth -a):
 Da la cobertura por cada base. El output son tres columnas: cromosoma, posición y cobertura.
Picard (CollectAlignmentSummaryMetrics): This tool takes a SAM/BAM file input and produces metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters. Nosotros tomamos el valor de PCT_Chimeras y lo incluímos en el genome results summary. 


Funciones que hemos programado (integrados en el script principal como funciones): 
parse_bam_header():
Lee el header del BAM y saca información básica como la versión del formato, el orden del sort, la plataforma usada, la muestra (sample) y las herramientas usadas durante el alineamiento.


calculate_mapq_median()
Se centra en calcular cómo varía la calidad de mapeo (MAPQ) a lo largo del genoma. Toma el cromosoma, la posición, y la calidad de mapeo de cada lectura, agrupa en ventanas de 3 millón de pares de bases, y calcula la mediana por ventana. Al final, guarda el resultado en un json.


MÓDULOS DE MULTIQC 
Tablas personalizadas
→ Genome Results Summary
MultiQC trae por defecto una tabla general de estadísticas (General Stats Table), pero vimos que no era capaz de parsear toda la información que nos interesa de los archivos generados por herramientas como Qualimap. Por eso decidimos eliminar esa tabla por defecto y crear la nuestra propia, usando un script que lee el archivo genome_results.txt y extrae las métricas que realmente necesitamos:
Variables que recogemos:
Sample
File name 
Número total de lecturas
Lecturas mapeadas
Lecturas mapeadas en par (ambas)
Lecturas mapeadas como singletons
% de lecturas duplicadas (flagged)
Tamaño medio del insert
Calidad de mapeo media
Porcentaje de GC
Tasa de error general (mismatches)
Cobertura media
Desviación estándar de la cobertura
PCT_chimeras (picard)



→ BAM header information
Los archivos BAM traen mucha información en su cabecera, y como esa parte también nos interesa, usamos un script para extraer los siguientes campos y mostrarlos como tabla en el reporte:
Sample name 
Versión del formato BAM
Orden del sort
Ensamblado usado
Plataforma de secuenciación
Herramientas utilizadas



Cómo crear un módulo nuevo en MultiQC
Objetivo:
 Crear un módulo propio en MultiQC que sea capaz de leer nuestros archivos bam_header_info.txt y genome_results.txt, y generar dos tablas nuevas que se incluyan en el informe final de MultiQC.
Pasos que seguimos:
Clonamos el repositorio oficial:

 bash
CopyEdit
git clone https://github.com/MultiQC/MultiQC.git

Nos movimos a la versión que necesitábamos:

 bash
CopyEdit
git checkout v1.25.2

Creamos un branch nuevo (para no tocar el código original):

 bash
CopyEdit
git checkout -b Aurora

Seguimos la documentación oficial para desarrollo de módulos:
 https://docs.seqera.io/multiqc/development/modules


Creamos la carpeta del nuevo módulo dentro de multiqc/modules y la llamamos tablemaker. La estructura queda así:


multiqc/
├── modules/
│   └── tablemaker/
│       ├── __init__.py
│       ├── tablemaker.py
│       └── tests/
│           └── test_tablemaker.py
Creamos archivos clave dentro de la carpeta que acabamientos de crear
→ __init__.py
En este archivo escribimos: 
from .tablemaker import MultiqcModule
__all__ = ["MultiqcModule"]

→ pyproject.toml
Ahí le decimos a MultiQC que nuestro módulo existe, añadiendo en la sección [project.entry-points."multiqc.modules.v1"]:
tablemaker = "multiqc.modules.tablemaker:MultiqcModule"
Luego reinstalamos MultiQC localmente con:
pip install -e .
Definimos el orden de aparición
Para que las tablas aparezcan primero en el reporte, editamos config_defaults.yaml y añadimos tablemaker al principio de module_order.
module_order:
  - tablemaker
  - qualimap
  - ...

Podemos añadir tests
Aunque MultiQC no requiere tests, vimos que algunos módulos los tienen, como:
samtools/tests


picard/tests


umitools/tests


Nos parece una buena idea crear nuestros propios tests en el futuro (TO DO), para asegurarnos de que todo sigue funcionando cuando haya cambios.
Definimos el tipo de módulo
MultiQC tiene “helpers” para crear distintos tipos de módulos. Nosotros estamos usando uno del tipo tabla. El código de base es este:
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
import logging

log = logging.getLogger(__name__)
log.info("Hello World!")

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Table Maker",
            anchor="tablemaker",
            info="Parse of information from genome_results.txt from Qualimap and BAM header"
        )
Cómo detecta nuestro módulo los archivos
En el archivo search_patterns.yaml definimos qué archivos necesita leer el módulo:
yaml
CopyEdit
tablemaker/header:
  fn: "bam_header_info.txt"
  num_lines: 25

tablemaker/genome:
  fn: "genome_results.txt"
  num_lines: 25

⚠️ Notamos que MultiQC ya intentaba parsear genome_results.txt desde Qualimap, y por eso nuestro módulo no lo encontraba. Para arreglarlo, quitamos ese archivo de la sección de qualimap en el mismo search_patterns.yaml. De paso, esto eliminó la General Stats Table por defecto, que era algo que queríamos hacer, así que, perfecto.

Parseo de archivos
Creamos una función para leer los archivos genome_results.txt y bam_header_info.txt
def parse_tabbed_key_value(file_contents):
    data = {}
    for line in file_contents.strip().splitlines():
        parts = line.strip().split("\t")
        if len(parts) == 2:
            key, value = parts
            data[key.strip()] = value.strip()
    return data

Esto nos permite parsear rápidamente archivos como el genome_results.txt o el bam_header_info.txt, y con esos datos creamos las tablas. 
Después, tenemos que empezar un nuevo archivo dentro de la carpeta de nuestro modulo (tablemaker en este caso) y crear un script llamado tablemaker.py donde finalmente se hace el proceso de parsear los archivos y crear el html. 

Cómo ejecutar
Una vez hecho todo lo anterior, generamos el informe con:
cd /home/aumoreno/Desktop/codigo_multiqc/MultiQC
multiqc /home/aumoreno/Desktop/bam/final_github/bam_analysis_results --force



2. MultiQC graphics from qualimap
Para empezar, entender qué es qualimap, es una herramienta que genera estadísticas generales sobre archivos bam. Esta herramienta genera diferentes archivos .txt con información relevante, entre ellos se encuentran (docu de este sitio web): 
Insert size histogram: Histogram of insert size distribution. To construct the histogram all collected insert size values (number of read alignments of a certain insert size) are seprated in bins. The default number of bins is 50. The number of reads in the bin will be the sum from insert size of the bins. The detailed values for insert are reported in raw data report.
Genome fraction coverage: Provides a visual way of knowing how much reference has been sequenced with at least a given coverage rate. This graph should be interpreted as in this example:
If one aims a coverage rate of at least 25X (x-axis), how much of reference (y-axis) will be considered?
Mapped reads clipping profile: Represents the percentage of clipped bases across the reads. The clipping is detected via SAM format CIGAR codes ‘H’ (hard clipping) and ‘S’ (soft clipping). In addition, the total number of clipped reads can be found in the report Summary. The plot is not shown if there are no clipped-reads are found. 
Mapping quality histogram: Histogram of the number of genomic locations having a given mapping quality. To construct the histogram mean mapping quality is computed at each genome position with non-zero coverage and collected. According to Specification of the SAM format the range for the mapping quality is [0-255].
Homopolymer indels: This bar plot shows separately the number of indels that are within a homopolymer of A’s, C’s, G’s or T’s together with the number of indels that are not within a homopolymer. Large numbers of homopolymer indels may indicate a problem in a sequencing process. An indel is considered homopolymeric if it is found within a homopolymer (defined as at least 5 equal consecutive bases). Owing to the fact that Qualimap works directly from BAM files (and not from reference genomes), we make use of the CIGAR code from the corresponding read for this task
Mapped reads nucleotide content: This plot shows the nucleotide content per position of the mapped reads.

A pesar de que el módulo de qualimap existe en multiqc nos hemos encontrado con muchos inconvenientes a la hora de convertir sus archivos a gráficos. Para empezar, no todos los gráficos que nos interesan están programados por multiqc para ser parseados, y además, no siempre lo hacen de la mejor forma. En versiones posteriores de la herramienta estos errores se mejoran, pero el gráfico mapped reads per nucleotide content no funciona en otras versiones. Por lo tanto, debemos llevar a cabo ciertas modificaciones. Se realizan en el script BAM_finalize_2.py, es requisito eliminar el encabezado de todos los archivos, y por otro lado, es requisito, añadirle al archivo homopolymer_indels labels llamadas Labels y Number of Indels, para el futuro gráfico. También es necesario añadir el sufijo _mqc.txt a los archivos para que multiqc los detecte. El único gráfico no apto para la interfaz de MultiQC es mapped reads nucleotide content, el cual necesita ser modificado completamente, en concreto, se manipula el título y se transponen las columnas.
Próximamente: Ahora mismo, lo que tenemos es un gráfico de tipo Line, y en el caso de los histogramas, no se observan barras, cuando se podrían observar, para ello, se puede forzar a que aparezcan simplemente añadiendo nombres a las columnas separados por una tabulacion, ejemplo
 Insert size 	number of reads 
Pero el gráfico se da la vuelta: 


y se observa de esta forma, lo comentaremos en proximas reuniones. 
Opciones de mejorar esto: 
Modificar el código del módulo de qualimap. 
Modificar el archivo multiqc_config.yaml para que incluya todos los datos presente en los archivos que queremos plotear y a partir de ahí configurar qué queremos que sea el eje x, el eje Y, los nombres de las labels… (custom_content)
Intentar modificar el contenido de los archivos para que se adapte a lo que el modulo custom content espera por defecto, y sin tener que modificar ni el archivo config, ni el modulo de multiqc para qualimap se parseen “solos”, entiendo que el orden de las columnas debe de modificarse, debe de añadirse una sección de configuración estilo yaml al principio del archivo… 
Tras comprobar esta información y testar (demasiadas horas) no se ha conseguido, y en la propia docu de la herramienta dejan claro de que no se pueden crear histogramas, así que nos quedamos cómo estamos (gráfico de tipo Line) 
3. Módulo: custom graphics across reference
Se querían realizar dos gráficos donde se pudiese visualizar la calidad de mapeo y la cobertura respectivamente a lo largo del genome secuenciado. Para ello, se necesita tener como referencia el cromosoma y la posición. 
COVERAGE ACROSS REFERENCE

Se necesitaba obtener la cobertura por posición y cromosoma de todo el genoma. Para ello se testaron dos herramientas:

Tools Compared: bedtools, samtools
Version(s): Samtools Version: 1.13 (using htslib 1.13+ds), bedtools Version:   v2.30.0
Date of Testing: April 2025
Environment: Linux Ubuntu 22.04

Datasets Used 

1000 genomes project: low coverage bam file, de la fase 3 del proyecto. 
Se puede descargar  aquí
En mi ordenador local se encuentra en la siguiente dirección:
/home/aumoreno/Desktop/bam/coverage/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 
BAM de cromosoma 1 que obtuve a partir de biobox para casos de emergencia como este. Lo obtuve seleccionando el cromosoma 1 a partir de un BAM de WGS, el BAM original exacto no lo recuerdo, pero lo obtuve de este dataset. En local lo tengo ubicado en /home/aumoreno/Desktop/bam/mediana/chr1.bam
Tool configuration: 
Commands:
SAMTOOLS:

samtools depth -a input.bam > depth.txt 

BEDTOOLS:

bedtools genomecov -d input.bam > depth.txt

Benchmarking and performance: 

Bedtools tarda una gran cantidad de tiempo, del orden de 3/4 horas para un bam de genoma completo en generar la cobertura, cuando samtools ha demostrado tardar unos 10 minutos, es por ello que he seleccionado samtools como herramienta para obtener esta información. Ambos generan archivos muy pesados, del orden de 40gb para un BAM de low coverage, por lo que sería importante hacer tests en bam muy pesados y de mucha información.

La herramienta samtools con el comando depth -a es la seleccionada, el output se manipula con un script (adaptado como función dentro del script BAM_finalize_2.py), que generará intervalos de 3M de pb, y finalmente dará el valor de la mediana para ese intervalo. El output será un archivo json que debe plotearse de forma que el archivo final tenga la siguiente forma: 
En lugar de guardar en disco la salida intermedia de herramientas como samtools o pysam en el clúster de Starlife, se ha optado por redireccionar dicha salida directamente al script que genera los intervalos genómicos de 3 millones de pares de bases (3Mb). Este script se encuentra en local bajo la ruta:
/home/aumoreno/Desktop/bam/coverage/intervalos.py
Características de la función:
Ventanas genómicas acumuladas:
 El genoma se divide en ventanas de 3Mb no por cromosoma individual, sino de forma acumulativa. Es decir, la posición inicial del cromosoma 2 se suma al tamaño del cromosoma 1. La posición final del cromosoma 3 equivale a la suma de los tamaños de los cromosomas 1, 2 y 3. Esto garantiza una representación lineal continua del genoma en el eje X de los gráficos.


Longitudes cromosómicas desde el header (LN):
 Cuando el archivo BAM contiene header, las longitudes de los cromosomas se extraen directamente desde las líneas @SQ, campo LN, si no tiene header, estas longitudes no se toman y el gráfico no funciona. 


Estandarización de nombres de cromosomas:
 Dado que los BAM pueden contener nombres de referencia no estándar o inconsistentes (por ejemplo, “1”, “MT”, “M”, “chr1”, etc.), se utiliza un archivo de alias (chrnames.json) que traduce estos nombres a un formato canónico tipo "chrN", "chrX", "chrM", etc, para que las labels de los cromosomas en el gráfico no cambien


Problema: muchos archivos (~700.000 BAM/CRAM, según los datos de Ángel, se entiende que se están teniendo en cuenta archivos no released) no contienen header, y por lo tanto no es posible obtener las longitudes de los cromosomas con .get_reference_length(). Tampoco hay garantía de que los nombres de referencia sean estándar.
Solución: Se ha decidido no modifcar el script a pesar de esto, en el futuro no se podrán realizar submission de BAM sin header, por lo que se entiende que el problema no se agravará.


Resultado final:
 El script produce un archivo .json donde cada línea contiene tres columnas: el cromosoma (formato chrN), la posición de inicio global de la ventana, y la mediana (de cobertura o MAPQ). Esto se puede usar para representar visualmente la distribución de calidad o cobertura sobre el eje genómico completo de forma continua, incluso si solo hay un subconjunto de cromosomas secuenciados.

MEDIANA ACROSS REFERENCE

Script se encuentra en /home/aumoreno/Desktop/bam/mediana/1000_genomes/mt.py

El siguiente gráfico es la mediana de la calidad a lo largo del genoma. No se encontró ninguna herramienta que permitiese el cálculo de la misma, intenté utilizar el output que genera qualimap en mapping_quality_across_reference.txt, pero al no obtener mucha información sobre cómo qualimap generaba esos datos, finalmente decidí llevar a cabo yo misma un script que tomase las columnas del BAM correspondientes al cromosoma, a la posición y a la calidad de mapeo en esa misma posición. Hay que tener en cuenta que en la mayoría de posiciones hay más de un read solapando en la misma posición, por lo que en primer lugar se calcula la mediana de todos los reads que coinciden en la misma posición del mismo cromosoma, y una vez obtenido esto, se calculan intervalos de 3M de pares de bases para todo el genoma, y la mediana correspondiente a cada intervalo, con el objetivo de representarlo. Se intentó imitar el mismo gráfico que genera qualimap, y se obtuvieron los siguientes resultados: 


Imagen generada a partir del script diseñado por mí. 





Imagen obtenida de ejecutar qualimap.

Como vemos, la gran mayoría de picos coinciden, me fijé en el pico del cromosoma 8 a 9, que llega a 0, y como en el gráfico generado por mí este pico no coincidía. Si observamos los archivos sobre los que se basan ambos gráficos, vemos que el pico coincide en la posición, el problema es que qualimap interpreta de forma errónea el final del cromosoma 8 y el inicio del cromosoma 9. Por lo que confirmamos que hemos conseguido ejecutar el gráfico de forma correcta. 

Ejemplo generado por mi: 



qualimap: 



El equipo de dev team ha desarrollado este mismo gráfico de forma interactiva, y así tanto la calidad como el coverage sea fácilmente observable a lo largo del cromosoma. 
¿Cómo lo han hecho?
Creando dos módulos, uno para cada tabla, llamados genomewide_mapq y genomewide_coverage. Cada uno plotea un json con los valores de calidad y cobertura respectivamente a lo largo de cromosomas. Se puede hacer zoom y scrolling. 

¿Cómo se ejecuta multiqc? 
Debes tener una carpeta, llamada (bam_analysis_results) con todos estos contenidos: 

bam_header_info.txt 
 coverage_median_3Mb_complete.json  genome_results.txt     
mapq_median_3Mb_complete.json      
 read_distribution.txt
css                                
multiqc_config.yaml          
 raw_data_qualimapReport

y dentro de raw_data_qualimapReport, los siguientes contenidos: 

genome_fraction_coverage_mqc.txt 
 homopolymer_indels_mqc.txt 
 insert_size_plot_mqc.txt  
mapped_reads_nucleotide_content_mqc.txt  
mapping_quality_plot_mqc.txt
Fuera de la carpeta bam_analysis_results ejecutar el siguiente comando. 
multiqc . -force -c bam_analysis_results/multiqc_config.yaml

Modificaciones de MultiQC a posteriori

Eliminar el botón de porcentajes del gráfico de homopolymer indels. 
Eliminar el label del número de samples del gráfico de homopolymer indels y mapped nucleotide content. 
Eliminar el label del día y la hora en la que se generó el report y el path del archivo. 
