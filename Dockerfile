FROM ubuntu:latest

# Update and upgrade packages
RUN apt-get update && apt-get -y upgrade

# Install Python and virtualenv
RUN apt-get install -y python3-pip python3-venv

# Create a directory for the user and set up a virtual environment
RUN mkdir -p /home/user
RUN python3 -m venv /home/user/venv

# Copy requirements.txt to the working directory
COPY requirements.txt /home/user/requirements.txt

# Install Python dependencies in the virtual environment
RUN /home/user/venv/bin/pip install -r /home/user/requirements.txt

# Conditionally copy and install NUPACK based on the architecture
RUN uname -m | grep -q aarch64 && \
    (COPY nupack-4.0.1.8/package/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_aarch64.manylinux2014_aarch64.whl /home/user && \
    /home/user/venv/bin/pip install /home/user/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_aarch64.manylinux2014_aarch64.whl) || \
    (uname -m | grep -q x86_64 && \
    (COPY nupack-4.0.1.8/package/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl /home/user && \
    /home/user/venv/bin/pip install /home/user/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl))

# Set the default command to use the virtual environment
CMD ["/bin/bash", "-c", "source /home/user/venv/bin/activate && exec /bin/bash"]
