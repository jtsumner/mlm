from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from os.path import join
FTP = FTPRemoteProvider(username="email", password="pswd", port=21)
url1 = "ftp.box.com/JackSumner_Data/First_Batch_of_Reads/"
url2 = "ftp.box.com/JackSumner_Data/Second_Batch_of_Reads/"

#### Rules to download files from box for FTP ####

rule ftp_first:
    input:
        FTP.remote(join(url1, "{file_f}"), immediate_close=True)
    output:
        temp("../data/reads/{file_f}")
    shell:
        "mv {input} {output}"

rule ftp_second:
    input:
        FTP.remote(join(url2, "{file_s}"), immediate_close=True)
    output:
        temp("../data/reads/{file_s}")
    shell:
        "mv {input} {output}"
