
sys="IEEE14"

lvls1=["1e0","1e1","1e2","1e3","1e4","m1e0","m1e1","m1e2","m1e3","m1e4"]
lvls2=["5","10","15","20","m5","m10","m15","m20"]


path="Teste3/"+sys+"/"


for lvl in lvls2:

    input_file = path+"conds_"+lvl+".csv"
    output_prefix = path+ "conds_"+lvl

    string=["A","B"]


    with open(input_file, "r") as f:
        file_num = 0
        output_file = open(output_prefix + string[file_num] + ".csv", "w")
        for line in f:
            if "Estimador 2" in line:
                output_file.close()
                file_num += 1
                output_file = open(output_prefix + string[file_num] + ".csv", "w")
            output_file.write(line)
        output_file.close()