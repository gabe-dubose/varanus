
#function to write database to file
def write_database_file(database, filename):
    import json

    with open(filename, 'w') as outfile:
        json.dump(database, outfile)

#function to write variant annotation to output file
