#!/usr/bin python
import csv
# Do the reading
file1 = open('2_gene_extended_model/data/model_codons.csv', 'rb')
reader = csv.reader(file1,delimiter='\t')
new_rows_list = []
for row in reader:
    print(row)
    new_row = [row[0], row[1]]
    #Decrement the value of the slow codon from 10->1->.4->.2 so it gets faster
    #if new_row[0] == 'TAC':
    #    if float(new_row[1])   == 0.111:
    #        new_row[1] = float(0.10)
    #    elif float(new_row[1]) == 0.133:
    #        new_row[1] = float(0.111)
    #    elif float(new_row[1]) == 0.20:
    #        new_row[1] = float(0.133)
    #    elif float(new_row[1]) == 0.40:
    #        new_row[1] = float(0.20)
    #    elif float(new_row[1]) == 1.00:
    #        new_row[1] = float(0.40)
    #    elif float(new_row[1]) == 10.0:
    #        new_row[1] = float(1.00)
    #    elif float(new_row[1]) == .10:
    #        new_row[1] = float(10.00)
    #    else:
    #        print(new_row[1])
    #        raise ValueError('Unspecified Codon Values')
    if new_row[0] == 'TAT': #Changed so that we are slowing down our fast codon
        if float(new_row[1])   == 0.90:
            new_row[1] = float(1.00)
        elif float(new_row[1]) == 0.75:
            new_row[1] = float(0.90)
        elif float(new_row[1]) == 0.50:
            new_row[1] = float(0.75)
        elif float(new_row[1]) == 0.25:
            new_row[1] = float(0.50)
        elif float(new_row[1]) == 0.10:
            new_row[1] = float(0.25)
        elif float(new_row[1]) == 0.01:
            new_row[1] = float(0.10)
        elif float(new_row[1]) == 1.00: #Reset
            new_row[1] = float(0.01)
        else:
            print(new_row[1])
            raise ValueError('Unspecified Codon Values')
    new_rows_list.append(new_row)
file1.close()   # <---IMPORTANT
# Do the writing
file2 = open('2_gene_extended_model/data/model_codons.csv', 'wb')
writer = csv.writer(file2,delimiter='\t')
writer.writerows(new_rows_list)
file2.close()
