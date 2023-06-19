import sys,os


def getReadEndpoint(start,mapinfo):
    currnum = ""
    totallen = 0
    for char in mapinfo:
        if char.isdigit():
            currnum += char
        else:
            currnum = int(currnum)
            if char == "D":
                totallen += currnum
            elif char == "I":
                pass
            elif char == "M":
                totallen += currnum
            elif char in ["S", "H"]:
                pass
            else:
                raise Exception
            currnum = ""
    return start+totallen-1

line = sys.stdin.readline()
readh = {}
headers = ["@HD", "@PG", "@RG", "@SQ"]
while line:
    if not line[:3] in headers:
        line = line.strip().split("\t")
        read,flag,c,pos1,mapqual,mapinfo,c2,pos2,isize,reads,quals = line[:11]
        mapqual = int(mapqual)
        isize = int(isize)
        intFlag = int(flag)
        flag = bin(intFlag)
        #the read is paired, primary, and passes qual checks and is not a PCR/optical duplicate
        if flag[-1] == "1" and intFlag <= 512:
            #make sure both reads are mapped
            if flag[-3] == "0" and flag[-4] == "0":
                pos1 = int(pos1)
                pos2 = int(pos2)
                if flag[-5] == "0":
                    strand = "+"
                else:
                    strand = "-"
                if flag[-6] == "0":
                    strand2 = "+"
                else:
                    strand2 = "-"
                if c2 == "=" and pos1 < pos2 and (strand2+strand == "--" or strand2+strand == "++"):
                    if not readh.has_key(read):
                        readh[read] = ["",""]
                    s = pos1
                    e = getReadEndpoint(s,mapinfo)
                    if readh[read][0] != "":
                        print read, flag, intFlag
                        raise Exception
                    readh[read][0] = (c,s,e,strand,mapqual,reads,quals)
                elif c2 == "=" and pos2 < pos1 and (strand2+strand == "--" or strand2+strand == "++"):
                    #same chromosome and this one is on the right
                    if not readh.has_key(read):
                        readh[read] = ["",""]
                    s = pos1
                    e = getReadEndpoint(s,mapinfo)
                    if readh[read][1] != "":
                        raise Exception
                    readh[read][1] = (c,s,e,strand,mapqual,reads,quals)
                if readh.has_key(read) and readh[read][0] != "" and readh[read][1] != "":
                    lcoords,rcoords = readh[read]
                    c1,ls,le,strand1,mapqual1,reads1,quals1 = lcoords
                    c2,rs,re,strand2,mapqual2,reads2,quals2 = rcoords
                    if c1 != c2:
                        raise Exception
                    if rs < ls:
                        print read,lcoords,rcoords
                        raise Exception
                    span = (rs - le) - 1
                    strands = strand1+strand2
                    if strands == "-+":
                        print "\t".join([str(x) for x in [read,c1,ls,le,strand1,mapqual1,reads1,quals1,rs,re,strand2,mapqual2,reads2,quals2]])
                        del readh[read]
    line = sys.stdin.readline()
