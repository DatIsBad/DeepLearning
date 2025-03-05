from e1_DNAtoProtein import DNAtoProtein as e1
from e2_ProteinToTable import ProteinToTable as e2
from e3_createDatabase import Database as e3
from e4_evaluate import evaluate as e4

db = e3.Database()


#create database + tables if they don't exit
db.create_database()

#reset database
db.reset_database()

#erase tables
#db.dropAll()

#get samples
sample_list = e2.openFile(e2.path1)

#print few samples
i = 0
for sample in sample_list:
    print(sample)
    i = i + 1
    if i > 10:
        break


db.insert_data(sample_list)