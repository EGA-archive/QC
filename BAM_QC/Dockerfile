# Lightweight Python image
FROM python:3.11.18-slim

# Install system-level dependencies
RUN apt-get update && apt-get install -y \
    samtools \
    openjdk-17-jre \
    git \
    wget \
    ca-certificates \
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
# Security fix for eMudhra/emSign CA incident
RUN python -m pip install \
    "pip @ https://files.pythonhosted.org/packages/3f/06/3e4420302363fbb71be3c7f35ce042f54bfca330e3e66b3a2e3c4cca579f/pip-26.2.2-py3-none-any.whl#sha256=fd7f4dfa8190e014013b9f08f3e973f4e4c26ffd179a204895533a7254813a70"

# Entry point that executes your wrapper
ENTRYPOINT ["python", "run/wrapper.py"]
