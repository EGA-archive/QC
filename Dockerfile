# Lightweight Python image
FROM python:3.10-slim

# Install system-level dependencies
RUN apt-get update && apt-get install -y \
    samtools \
    openjdk-17-jre \
    git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy the entire project into the container
COPY . /app

# Ensure the Qualimap binary and scripts are executable (already present in run/ and output/)
RUN chmod +x run/qualimap
RUN chmod +x run/BAM_pipeline_2.py output/BAM_finalize_2.py

# Install Python dependencies (including your custom MultiQC)
RUN pip install --no-cache-dir -r requirements.txt

# Entry point that executes your wrapper
ENTRYPOINT ["python", "run/wrapper.py"]
