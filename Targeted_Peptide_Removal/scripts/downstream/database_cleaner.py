import os

def delete_temp_database_files():
    # Define the directory
    database_dir = "database"

    # Check if the directory exists
    if not os.path.exists(database_dir):
        print(f"Directory '{database_dir}' does not exist.")
        return

    # Get list of temp_database files
    files_to_delete = [f for f in os.listdir(database_dir) if f.startswith("temp_database")]

    # Delete the files
    for file_name in files_to_delete:
        file_path = os.path.join(database_dir, file_name)
        try:
            os.remove(file_path)
            print(f"Deleted: {file_name}")
        except Exception as e:
            print(f"Failed to delete {file_name}: {e}")

# Execute the function
delete_temp_database_files()
