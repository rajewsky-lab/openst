# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Update and install necessary packages
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3-pip \
    python3-venv \
    libxcb-xinerama0 \
    libnss3 \
    libxcomposite1 \
    libxcursor1 \
    libxi6 \
    libxtst6 \
    libxrandr2 \
    libasound2 \
    libatk1.0-0 \
    libatk-bridge2.0-0 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgtk-3-0

# Install openst
RUN pip3 install openst

# Set up X11 for GUI applications
ENV DISPLAY=:0

# Create a directory for the application
WORKDIR /app

# Start the terminal when the container runs
CMD ["/bin/sh", "-c", "bash"]