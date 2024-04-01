list = []
with open('..\data\Protein_f\sequence.txt','r',encoding='utf-8') as fr:
    for line in fr:
        if line == '\n' : continue
        list.append(line)
        if line.__contains__('VERSION') :
            fileName = '..\data\Protein_f\phage-gb\\' + line[12:-3] + '.txt'
        if line.__contains__('//') and len(line) <= 4 :
            with open(fileName,'w',encoding='utf-8') as fw:
                fw.writelines(list)
            list = []

