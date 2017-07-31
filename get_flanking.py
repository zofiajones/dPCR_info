import pandas as pd
import tabix
import pysam

file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

file1=open('flanking_snps.txt','r')
lines1=file1.readlines()

file2=open('flanking_regions.txt','w')

base0=''
for line in lines1:
	m=line.strip().split()
	print(m)
	if len(m)==2:
		if base0:
			#print(middle,middle_end,origin,line0)
			#print(base0)
			vars=sorted(vars)
			ind=vars.index(origin)
			n=base0.split('[')
			print(ind,len(vars),len(n))
			p1=''
			p1=p1+n[0]+'['
			for item in range(1,ind+1):
				p1=p1+n[item]+'['
			p2=''
			p2=p2+n[ind+1].split(']')[1]+'['
			for item in range(ind+2,len(n)):
				print(item)
				p2=p2+n[item].split(']')[0]+']'
			if (ind+1) != (len(n)-1):
			#if ']' in n[len(n)-1]:
				p2=p2+n[len(n)-1].split(']')[1]+']'
			file2.write(line0+'\n')
#			file2.write(base0+'\n')
			file2.write(p1[:-1]+'\n')
			file2.write('['+n[ind+1].split(']')[0]+']'+'\n')
			file2.write(p2[:-1]+'\n\n')
			#file2.write(base0[middle:middle_end+1]+'\n')
			#file2.write(base0[middle_end+1:]+'\n')
			#file2.write(str(found_origin)+'\torigin\n')
		loc=m[0]+':'+ str(int(m[1])-150)+'-'+str(int(m[1])+150)
		base0=file.fetch(region=loc)
		#print(len(base0))
		origin=int(m[1])
		vars=[]
		middle=150
		middle_end=0
		found_origin=0
	if len(m)==4:
		A=list(set(m[2]))
		B=list(set(m[3]))
		ref=m[2]
		alt=m[3]
		if int(m[1])==origin:
			line0=line.strip()
		if A[0]=='O':
			m[1]=str(int(m[1])-1)
			diff=int(m[1])-origin+middle
			p1=base0[:diff+1]
			p2=base0[diff+1:]
			print(base0[middle])
			print('base0',base0[diff-5:diff+2],base0[diff+1:diff+6])
			base0=p1+'[*/'+alt+']'+p2
#			print(base0)
#			print(len(p1),len(p2))
			if (int(m[1])+1) < origin:
				add='[*/'+alt+']'
				middle=middle+len(add)-1
			vars.append(int(m[1])+1)
		elif B[0]=='O':
			m[1]=str(int(m[1])-1)
			diff=int(m[1])-origin+middle
			print('base0',base0[diff-5:diff+2],base0[diff+2:diff+10])
			p1=base0[:diff+1]
			p2=base0[diff+len(ref)+1:]
			base0=p1+'['+ref+'/*]'+p2
			#print(len)
#			print(len(p1),len(p2))
			if (int(m[1])+1) < origin:
				add='['+ref+'/*]'
				print('len',len(add))
				middle=middle+len(add)-1
				diff=int(m[1])-origin+middle
				print('diff',diff)
				print(base0[middle])
			vars.append(int(m[1])+1)
		else:
			diff=int(m[1])-origin+middle
			print('origin',base0[diff])
			print('origin',base0[diff-5:diff],base0[diff+1:diff+5])
			p1=base0[:diff]
			p2=base0[diff+1:]
			base0=p1+'['+ref+'/'+alt+']'+p2
#			print(base0)
#			print(len(p1),len(p2))
			if int(m[1]) < origin:
				add='['+ref+'/'+alt+']'
				middle=middle+len(add)-1
				print(base0[middle])
			vars.append(int(m[1]))
		loc=m[0]+':'+ str(int(m[1]))+'-'+str(int(m[1]))
		base1=file.fetch(region=loc)
		if m[2] == base1:
			print('check')
		else:
			print('not check',list(set(m[2])))

vars=sorted(vars)
ind=vars.index(origin)
n=base0.split('[')
print(ind,len(vars),len(n))
p1=''
p1=p1+n[0]+'['
for item in range(1,ind+1):
	p1=p1+n[item]+'['
p2=''
p2=p2+n[ind+1].split(']')[1]+'['
for item in range(ind+2,len(n)):
	print(item)
	p2=p2+n[item].split(']')[0]+']'
if (ind+1) != (len(n)-1):
#if ']' in n[len(n)-1]:
	p2=p2+n[len(n)-1].split(']')[1]+']'
file2.write(line0+'\n')
#file2.write(base0+'\n')
file2.write(p1[:-1]+'\n')
file2.write('['+n[ind+1].split(']')[0]+']'+'\n')
file2.write(p2[:-1]+'\n\n')


file1.close()
file2.close()
