# Use the official Python base image
FROM python:3.9

# Set the working directory inside the container
WORKDIR /app

# Install libxrender1 package
RUN apt-get update && apt-get install -y libxrender1

# Copy the requirements file to the container
COPY requirements.txt .

# Install the dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application code to the container
COPY . .

# Set the command to run the Streamlit app
#CMD ["streamlit", "run", "main.py", "--server.enableCORS", "false", "--server.port", "8501", "--server.address", "0.0.0.0"]
CMD ["watchmedo", "auto-restart", "--directory", ".", "--pattern", "*.py", "--recursive", "--", "streamlit", "run", "main.py", "--server.port", "8502", "--server.address", "0.0.0.0"]
