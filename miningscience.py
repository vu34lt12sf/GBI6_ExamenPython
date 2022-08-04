from Bio import Entrez
import re
import pandas as pd
import matplotlib.pyplot as plt
import csv as csv
    


def download_pubmed(keyword): 
    """La función download:pubmed busca articulos con ID con las pabras claves, dentro de la base de datos pubmed"""
    
    from Bio import GenBank
    from Bio import Entrez
    from Bio import SeqIO
   
       
    Entrez.email = "gualapuro.moises@gmail.com"
    handle = Entrez.esearch(db="pubmed", 
                            term="Ecuador proteomics[Title/Abstract]",
                            usehistory="y")
    record = Entrez.read(handle)
    # generate a Python list with all Pubmed IDs of articles about Dengue Network
    id_list = record["IdList"]
    record["Count"]
    
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    handle = Entrez.efetch(db="pubmed",
                           rettype="medline", 
                           retmode="text", 
                           retstart=0,
    retmax=543, webenv=webenv, query_key=query_key)
    
    out_handle = open("data/ECmetagenome_pubs.txt", "w")
    data = handle.read()
    handle.close()
    out_handle.write(data)
    out_handle.close()
    
    return

 
def mapscience(archivo):
    
    """La función map_science permite los paises que publican
    ese alrticuo mediante una grafica"""
    
    contentssc = re.sub(r'\s[\w._%+-]+@[\w.-]+\.[a-zA-Z]{1,4}','',archivo)
    contentsna = re.sub(r'\..\d.\,',',',contentssc)
    contentsnu = re.sub(r'\..\d.','',contentsna)
    x=contentsnu[1:].split('PMID-')
    
    Country1=[]
    for PMID in x:
        q=PMID.split('\n')
        for fila in q:
            w=fila.split(' ')
            if w[0] == 'AD':
                e=fila.split(',')
                Country1.append(e[-1])
    
    a=0
    Country2 =[0]*len(Country1)
    for lis in Country1:
        bytes(lis,encoding="utf8")
        if lis != '':
            w=lis
            if w[0] == ' ':
                w = re.sub (r'^\s','',w)
            if w[-1] == '.':
                w = re.sub (r'\.$','',w)
            w = re.sub (r'\.$','',w)
            w = re.sub (r'\s$','',w)
        Country2[a]=w
        a=a+1
    Countries=['Andorra','United Arab Emirates ','Afghanistan','Antigua and Barbuda','Anguilla','Albania','Armenia','Netherlands Antilles','Angola','Antarctica','Argentina','American Samoa','Austria','Australia','Aruba','Azerbaijan','Bosnia and Herzegovina','Barbados','Bangladesh','Belgium','Burkina Faso','Bulgaria','Bahrain', 'Burundi','Benin','Bermuda','Brunei','Bolivia', 'Brazil','Bahamas','Bhutan','Bouvet Island','Botswana','Belarus','Belize','Canada','Cocos [Keeling] Islands','Congo [DRC]','Central African Republic','Congo [Republic]', 'Switzerland',"Côte d'Ivoire",'Cook Islands','Chile','Cameroon','China','Colombia','Costa Rica','Cuba', 'Cape Verde','Christmas Island','Cyprus','Czech Republic','Germany','Djibouti','Denmark','Dominica','Dominican Republic','Algeria','Ecuador' ,'Estonia','Egypt','Western Sahara','Eritrea','Spain','Ethiopia','Finland','Fiji','Falkland Islands [Islas Malvinas]','Micronesia','Faroe Islands','France','Gabon', 'United Kingdom','Grenada','Georgia','French Guiana','Guernsey','Ghana','Gibraltar','Greenland','Gambia', 'Guinea','Guadeloupe','Equatorial Guinea','Greece','South Georgia and the South Sandwich Islands','Guatemala','Guam','Guinea-Bissau','Guyana','Gaza Strip','Hong Kong','Heard Island and McDonald Islands','Honduras','Croatia', 'Haiti','Hungary','Indonesia','Ireland' ,'Israel','Isle of Man','India','British Indian Ocean Territory','Iraq', 'Iran','Iceland','Italy','Jersey','Jamaica','Jordan', 'Japan','Kenya','Kyrgyzstan','Cambodia','Kiribati','Comoros','Saint Kitts and Nevis','North Korea','South Korea','Kuwait','Cayman Islands','Kazakhstan','Laos','Lebanon','Saint Lucia','Liechtenstein','Sri Lanka','Liberia','Lesotho','Lithuania','Luxembourg','Latvia' ,'Libya','Morocco','Monaco','Moldova','Montenegro','Madagascar','Marshall Islands','Macedonia [FYROM]','Mali','Myanmar [Burma]','Mongolia' ,'Macau','Northern Mariana Islands','Martinique','Mauritania','Montserrat','Malta','Mauritius','Maldives','Malawi','Mexico','Malaysia' ,'Mozambique','Namibia','New Caledonia','Niger','Norfolk Island','Nigeria','Nicaragua','The Netherlands','Norway','Nepal','Nauru', 'Niue','New Zealand','Oman','Panama','Peru','French Polynesia', 'Papua New Guinea','Philippines','Pakistan','Poland','Saint Pierre and Miquelon' ,'Pitcairn Islands','Puerto Rico','Palestinian Territories','Portugal','Palau','Paraguay','Qatar','Réunion','Romania', 'Serbia','Russia' ,'Rwanda','Saudi Arabia','Solomon Islands','Seychelles','Sudan','Sweden','Singapore','Saint Helena','Slovenia', 'Svalbard and Jan Mayen','Slovakia','Sierra Leone','San Marino','Senegal','Somalia','Suriname','São Tomé and Príncipe','El Salvador','Syria', 'Swaziland' ,'Turks and Caicos Islands','Chad','French Southern Territories','Togo','Thailand','Tajikistan','Tokelau','Timor-Leste','Turkmenistan' ,'Tunisia','Tonga','Turkey','Trinidad and Tobago','Tuvalu','Taiwan','Tanzania','Ukraine','Uganda','U.S. Minor Outlying Islands','United States of America','Uruguay','Uzbekistan','Vatican City','Saint Vincent and the Grenadines','Venezuela', 'British Virgin Islands','U.S. Virgin Islands','Vietnam','Vanuatu','Wallis and Futuna','Samoa','Kosovo','Yemen','Mayotte','South Africa','Zambia','Zimbabwe']
    Country3=Country2
    h=Countries
    f=len(h)
    Countryi=[0]*f
    k=0
    for elem in h:
        d=0
        for comp in Country3:
            if elem == str(comp):
                d=d+1
        Countryi[k]=d
        k=k+1
    
    Country4=[]
    Numb=[]
    o=0
    for line in Countryi:
        if str(line) != '0':
            Country4.append(line)
            m=Countries[o]
            Numb.append(m)
        o=o+1
        
    zip_coordinates = {}
    with open('Data/countries.txt') as f:
        csvr = csv.DictReader(f)
        for row in csvr:
            zip_coordinates[row['name']] = [float(row['latitude']),
                                            float(row['longitude'])]
    code = []
    long = []
    lat = []
    count = Country4
    for inte in Numb:
        if inte in zip_coordinates.keys():
            code.append(inte)
            lat.append(zip_coordinates[inte][0])
            long.append(zip_coordinates[inte][1])
    plt.scatter(long, lat, s = count, c= count)
    plt.colorbar()
    ard = dict(arrowstyle="->")
    plt.annotate('Singapore', xy = (103.819836, 1.352083), 
                 xytext = (110, 8), arrowprops = ard)
    plt.annotate('Guatemala', xy = (-90.230759, 15.783471), 
                 xytext = (-110, 20.4292), arrowprops= ard)
    plt.annotate('Belgium', xy = (4.469936, 50.503887),
                 xytext = (0, 55.25), arrowprops= ard)
    plt.annotate('United States of America', xy = (-95.712891,37.09024),
                 xytext = (-85.1106, 45.3736), arrowprops= ard)
    plt.annotate('Chile', xy = (-71.542969, -35.675147),
                 xytext = (-85.1106, -30.3736), arrowprops= ard)
    params = plt.gcf()
    plSize = params.get_size_inches()
    params.set_size_inches( (plSize[0] * 3, plSize[1] * 3) )
    plt.xlabel('Longitud')
    plt.ylabel('Latitud')
    plt.title('Map of Science')
    Mapa = plt.show()
    plt.savefig('img/Gráfica_primera búsqueda.jpg', dpi=300)
        
    return Mapa
   