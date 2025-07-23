# Lightweight Python image
FROM python:3.10-slim

# Install system-level dependencies
RUN apt-get update && apt-get install -y \
    samtools \
    openjdk-17-jre \
    git \
    wget \
    && apt-get clean && rm -rf /var/lib/apt/lists/*


# Set working directory
WORKDIR /app

# Copy the entire project into the container
COPY . /app

# Download picard.jar into the run/ directory
RUN wget -O run/picard.jar https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar

# Make scripts executable
RUN chmod +x run/qualimap_v2.3/qualimap \
             run/BAM_pipeline_2.py \
             output/BAM_finalize_2.py \
             run/wrapper.py

# Install Python dependencies (including your custom MultiQC)
RUN pip install --no-cache-dir -r requirements.txt

# Entry point that executes your wrapper
ENTRYPOINT ["python", "run/wrapper.py"]
