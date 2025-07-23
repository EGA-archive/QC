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

# Ensure the Qualimap binary is executable (already present in run/)
RUN chmod +x run/qualimap

# Install Python dependencies (including your custom MultiQC)
RUN pip install --no-cache-dir -r requirements.txt

# Entry point that executes your wrapper
ENTRYPOINT ["python", "run/wrapper.py"]
