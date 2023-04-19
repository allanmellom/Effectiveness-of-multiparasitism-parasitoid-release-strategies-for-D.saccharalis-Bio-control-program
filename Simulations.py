from random import SystemRandom 
import pickle as pck 
import math
import os
random=SystemRandom()
######
#P - Q - N - H or Q - H
#P- specialist parasitoid 
#Q- generalist parasitoid
#N- host of both parasitoid
#H- host of Q
######
#a1 - search efficiency of N by P (P -> N)
#a2 - search efficiency of N by Q (Q -> N)
#a3 - search efficiency of H by Q (Q -> H)
#lambda1 - population growth rate of N
#lambda2 - population growth rate of H
#H0=number of hosts needed to prevent female parasitoid to leave the patch ([PQ])
#h0=tolerance of others host ([NH])
#f0=tolerance of others female parasitoid ([PQ]) 
#taxa_disp -> dispersion rate of each population (PQNH)
#fracao_indv_migrante -> max dispersion rate of each population (PQNH/ 0Q0H)
#tempo_final -> last generation to simulate
#lin -> number of row
#col -> number of col
#pop_iniciais -> uma lista contendo as pop iniciais de cada ponto q tem pop inicial [[P,Q,N,H],[P,Q,N,H]], sendo q
	#pop_iniciais[0] vai ser estar no patch_iniciais[0]
#patch_iniciais -> uma lista contendo os patchs que ira iniciar com alguma pop [[l,c],[l,c]]
#viz-> é a lista dos vizinhos e distancias criada com a funcao criando_lista_viz(). É a primeira lista criada pela funcao citada.
#lista_raios-> list of radios in the choosen grid. It is the second list created by the function criando_lista_viz()
#TS= "total available searching time" ([P,Q1,Q2])
#TH = handling time ([P,Q1,Q2]) (tem dois Q pois o Q pega 2 host, então pode ter tempo de manuseio diferente para ambos. O q1 é para o host compartilhado e o Q2 para o host exclusivo do generalista.)
#porcentagem_popinihost= percentage of how many patches each host occupy at the start of simulations ([N,H])
#inicial_h = initial population of H (3)
#inicial_n = initial population of N (2)
#Pnumerico = how many adult especialist parasitoid emerge from a single main host
#Q1numerico = how many adult generalist parasitoid emerge from a single main host
#Q2numerico = how many adult generalist parasitoid emerge from a single alternative host



def modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,patches_sorteados,bOUm,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,borda_meio,loop,lib_extras,Pnumerico,Q1numerico,Q2numerico):
	qnd_onde_liberou=[]
	ocupacao_total=[[],[],[],[]] #criando a matriz de ocupacao ocupacao_total[pop][t]
	ocupacao_total_tempo=[0,0,0,0] #ocupacao_total_tempo[pop]

	media_regional=[[],[],[],[]] #criando a lista que vai acomodar a media regional media_regional[populacao][t]
	media_regional_tempo=[0,0,0,0] #criando a lista que vai acomodar a media regional media_regional[populacao][t]

	grid=lin*col #criando o valor de grid para poder fazer a media regional

	pop_migrante_tempo=[0,0,0,0]
	media_regional_migracao=[[0],[0],[0],[0]]

	#criando os vetores dos dados
	vetores=criador_vetor(lin,col,tempo_final)
	g_all=vetores[0]
	pop_migrante=vetores[1]
	visitacao_patchs=vetores[2]
	imigracao_patchs=vetores[3]
	emigracao_patchs=vetores[4]

	###########
	#INICIO DINAMICA POPULACIONAL
	###########

	for t in range(tempo_final):
		for l in range(lin):
			for c in range(col):
				
				if t>0:
					tempo_dinamica=dinamica_iterativa(g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante,c,l,t,a1,a2,a3,lambda1,lambda2,TS,TH,H0,h0,f0,ocupacao_total,Pnumerico,Q1numerico,Q2numerico)
					#############
					#migração
					#############
					tempo2=criando_popmigrante_densidade_patch(l,c,t,H0,h0,f0,fracao_indv_migrante,g_all,pop_migrante,emigracao_patchs)#criei a pop migrante do tempo
					pop_migrante=tempo2[0]
					emigracao_patchs=tempo2[1]
					for i,j in enumerate(pop_migrante[l][c]): #estou diminuindo a pop do patch pela pop q vai migrar
						g_all[l][c][i][t]-=j[t]
				else:
					tempo_dinamica=dinamica_tempo0(g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante,c,l,t,patch_iniciais)
				g_all=tempo_dinamica[0]
				pop_iniciais=tempo_dinamica[1]
				ocupacao_total_tempo=tempo_dinamica[2]
				media_regional_tempo=tempo_dinamica[3]
				pop_migrante=tempo_dinamica[4]

		#saiu do loop pela matriz, agora, no mesmo tempo, vem o segundo loop pela matriz para poder migrar:
		if t!=0: #isso eh para nao ter migracao no tempo zero
			#loop de lin e col
			for l in range(lin):
				for c in range(col):
					tempo=sorteio_vizinhos(g_all,pop_migrante,taxa_disp,viz,l,c,t,visitacao_patchs,imigracao_patchs,lista_distancia)
					g_all=tempo[0]
					visitacao_patchs=tempo[1]
					imigracao_patchs=tempo[2]
					pop_migrante_tempo=[0,0,0,0]
					for i,j in enumerate(pop_migrante[l][c]): #i é posição, j elemento
						pop_migrante_tempo[i]+=pop_migrante[l][c][i][t]

			for i,j in enumerate(pop_migrante[l][c]): #pondo a media regional
				media_regional_migracao[i].append(pop_migrante_tempo[i]/grid)

		###Bloco de checar se as pop possuem recurso
		if t!=0: #para soh acontecer dps do tempo 0
			for l in range(lin):
				for c in range(col):
					g_all=conferir_recurso_parasitoides(g_all,l,c,t)
		
		

		
		
		if t!=0: #aqui atualiza a media regional e ocupacao
			for l in range(lin):
				for c in range(col):
					if g_all[l][c][0][t]>0:
						ocupacao_total_tempo[0]+=1
						media_regional_tempo[0]+=g_all[l][c][0][t]
					if g_all[l][c][1][t]>0:
						ocupacao_total_tempo[1]+=1
						media_regional_tempo[1]+=g_all[l][c][1][t]
					if g_all[l][c][2][t]>0:
						ocupacao_total_tempo[2]+=1
						media_regional_tempo[2]+=g_all[l][c][2][t]
					if g_all[l][c][3][t]>0:
						ocupacao_total_tempo[3]+=1
						media_regional_tempo[3]+=g_all[l][c][3][t]


		ocupacao_total[0].append(ocupacao_total_tempo[0])
		ocupacao_total[1].append(ocupacao_total_tempo[1])
		ocupacao_total[2].append(ocupacao_total_tempo[2])
		ocupacao_total[3].append(ocupacao_total_tempo[3])

		media_regional[0].append(media_regional_tempo[0]/grid)
		media_regional[1].append(media_regional_tempo[1]/grid)
		media_regional[2].append(media_regional_tempo[2]/grid)
		media_regional[3].append(media_regional_tempo[3]/grid)


		media_regional_tempo=[0,0,0,0]
		ocupacao_total_tempo=[0,0,0,0]
		
		if lib_extras==True: #aqui faz novas liberações só se tiver falado para fazer novas liberações
			if t!=0 and t!=(tempo_final-1): #aqui ta verificando se vai ter mais liberação
				tempo_BOL=False
				tempo=checando_necessidade_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,lista_patches,lista_borda,qnd_onde_liberou)
				g_all=tempo[0]
				qnd_onde_liberou=tempo[1]

				if tempo[2]==True:
					tempo_BOL=True


			if t!=0 and tempo_BOL==True: #ISSO REPETE PQ ELE PODE TER ADICIONADO MAIS BIXOS
				for l in range(lin):
					for c in range(col):
						if g_all[l][c][0][t]>0:
							ocupacao_total_tempo[0]+=1
							media_regional_tempo[0]+=g_all[l][c][0][t]
						if g_all[l][c][1][t]>0:
							ocupacao_total_tempo[1]+=1
							media_regional_tempo[1]+=g_all[l][c][1][t]
						#if g_all[l][c][2][t]>0:
						#	ocupacao_total_tempo[2]+=1
						#	media_regional_tempo[2]+=g_all[l][c][2][t]
						if g_all[l][c][3][t]>0:
							ocupacao_total_tempo[3]+=1
							media_regional_tempo[3]+=g_all[l][c][3][t]
				ocupacao_total[0][t]=ocupacao_total_tempo[0]
				ocupacao_total[1][t]=ocupacao_total_tempo[1]
				ocupacao_total[2][t]=ocupacao_total_tempo[2]
				ocupacao_total[3][t]=ocupacao_total_tempo[3]

				media_regional[0][t]=media_regional_tempo[0]/grid
				media_regional[1][t]=media_regional_tempo[1]/grid
				media_regional[2][t]=media_regional_tempo[2]/grid
				media_regional[3][t]=media_regional_tempo[3]/grid


				media_regional_tempo=[0,0,0,0]
				ocupacao_total_tempo=[0,0,0,0]



	
	salva_arquivos(zzzz,g_all,pop_migrante,ocupacao_total,media_regional,media_regional_migracao,visitacao_patchs,imigracao_patchs,emigracao_patchs,a1,a2,a3,lambda1,lambda2,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,qnd_onde_liberou)
	return
def criando_lista_patches_cada_hectare():
	lin=50
	col=50
	hectare1=[]
	hectare2=[]
	hectare3=[]
	hectare4=[]
	hectare5=[]
	hectare6=[]
	hectare7=[]
	hectare8=[]
	hectare9=[]
	hectare10=[]
	for l in range(lin):
		for c in range(col):
			entrou_em_qnts=0
			if l<10 and c<25:
				hectare1.append([l,c])
				entrou_em_qnts+=1
			elif l>9 and l<20 and c<25:
				hectare2.append([l,c])
				entrou_em_qnts+=1
			elif l>19 and l<30 and c<25:
				hectare3.append([l,c])
				entrou_em_qnts+=1
			elif l>29 and l<40 and c<25:
				hectare4.append([l,c])
				entrou_em_qnts+=1
			elif l>39 and l<50 and c<25:
				hectare5.append([l,c])
				entrou_em_qnts+=1
			elif l<10 and c>24:
				hectare6.append([l,c])
				entrou_em_qnts+=1
			elif l>9 and l<20 and c>24:
				hectare7.append([l,c])
				entrou_em_qnts+=1
			elif l>19 and l<30 and c>24:
				hectare8.append([l,c])
				entrou_em_qnts+=1
			elif l>29 and l<40 and c>24:
				hectare9.append([l,c])
				entrou_em_qnts+=1
			elif l>39 and l<50 and c>24:
				hectare10.append([l,c])
				entrou_em_qnts+=1
	arquivoviz=f'50x50_patches_por_hectare.txt'
	a=open(arquivoviz,"ab")		 
	dados_em_pck=pck.dumps([hectare1,hectare2,hectare3,hectare4,hectare5,hectare6,hectare7,hectare8,hectare9,hectare10])   
	a.write(dados_em_pck)
	a.close()
	return [hectare1,hectare2,hectare3,hectare4,hectare5,hectare6,hectare7,hectare8,hectare9,hectare10]
def criando_lista_patches_cada_borda():
	lin=50
	col=50
	hectare1=[]
	hectare2=[]
	hectare3=[]
	hectare4=[]
	hectare5=[]
	hectare6=[]
	hectare7=[]
	hectare8=[]
	hectare9=[]
	hectare10=[]
	for l in range(lin):
		for c in range(col):
			entrou_em_qnts=0
			if (l<10 and c==24) or (l<10 and c==0) or (l==0 and c<25) or (l==9 and c<25):
				hectare1.append([l,c])
				entrou_em_qnts+=1
			elif (l>9 and l<20 and c==24) or (l>9 and l<20 and c==0) or (l==10 and c<25) or (l==19 and c<25):
				hectare2.append([l,c])
				entrou_em_qnts+=1
			elif (l>19 and l<30 and c==24) or (l>19 and l<30 and c==0) or (l==20 and c<25) or (l==29 and c<25):
				hectare3.append([l,c])
				entrou_em_qnts+=1
			elif (l>29 and l<40 and c==24) or (l>29 and l<40 and c==0) or (l==30 and c<25) or (l==39 and c<25):
				hectare4.append([l,c])
				entrou_em_qnts+=1
			elif (l>39 and l<50 and c==24) or (l>39 and l<50 and c==0) or (l==40 and c<25) or (l==49 and c<25):
				hectare5.append([l,c])
				entrou_em_qnts+=1
			elif (l<10 and c==25) or (l<10 and c==49) or (l==0 and c>24) or (l==9 and c>24):
				hectare6.append([l,c])
				entrou_em_qnts+=1
			elif (l>9 and l<20 and c==25) or (l>9 and l<20 and c==49) or (l==10 and c>24) or (l==19 and c>24):
				hectare7.append([l,c])
				entrou_em_qnts+=1
			elif (l>19 and l<30 and c==25) or (l>19 and l<30 and c==49) or (l==20 and c>24) or (l==29 and c>24):
				hectare8.append([l,c])
				entrou_em_qnts+=1
			elif (l>29 and l<40 and c==25) or (l>29 and l<40 and c==49) or (l==30 and c>24) or (l==39 and c>24):
				hectare9.append([l,c])
				entrou_em_qnts+=1
			elif (l>39 and l<50 and c==25) or (l>39 and l<50 and c==49) or (l==40 and c>24) or (l==49 and c>24):
				hectare10.append([l,c])
				entrou_em_qnts+=1
	arquivoviz=f'50x50_patches_por_hectare_bordas.txt'
	a=open(arquivoviz,"ab")		 
	dados_em_pck=pck.dumps([hectare1,hectare2,hectare3,hectare4,hectare5,hectare6,hectare7,hectare8,hectare9,hectare10])   
	a.write(dados_em_pck)
	a.close()
	return [hectare1,hectare2,hectare3,hectare4,hectare5,hectare6,hectare7,hectare8,hectare9,hectare10]
def criando_lista_distancia(lin,col,distancia_maxima):
	lista_distancia=[]
	for l in range(lin):
		for c in range(col):
			if l==0 and c==0:
				continue
			else:
				distancia=(((0-l)**2)+((0-c)**2))**(1/2)
				if distancia<=distancia_maxima:
					lista_distancia.append(distancia)
	lista_distancia=sorted(set(lista_distancia)) #qnd o loop termina, eu retiro os valores repetidos e ordeno a lista
	return lista_distancia
def checando_existencia_lista_distancias(lin,col,distancia_maxima):
	arquivodist=f'{lin}x{col}_max{distancia_maxima}_lista_distancia.txt'
	try:
		distancia_file=open(arquivodist,"rb")
		distancia_pck=distancia_file.read()
		distancia_file.close()
		lista_distancia=pck.loads(distancia_pck)
	except FileNotFoundError:
		print("Não achou a lista de distância. Iniciando processo de criação.")
		lista_distancia=criando_lista_distancia(lin,col,distancia_maxima)
		a=open(arquivodist,"ab")		 
		dados_em_pck=pck.dumps(lista_distancia)   
		a.write(dados_em_pck)
		a.close()
	return lista_distancia
def checando_existencia_lista_viz_e_lista_de_patches(lin,col,distancia_maxima):
	arquivoviz=f'{lin}x{col}vizinhos.txt'
	try:
		viz_file=open(arquivoviz,"rb")
		viz_pck=viz_file.read()
		viz_file.close()
		viz=pck.loads(viz_pck)
	except FileNotFoundError:
		print("Não achou a lista de vizinhos. Iniciando processo de criação.")
		viz=criando_lista_viz(lin,col,distancia_maxima)
		a=open(arquivoviz,"ab")		 
		dados_em_pck=pck.dumps(viz)   
		a.write(dados_em_pck)
		a.close()
	######aqui vai checar borda
	arquivoviz=f'50x50_patches_por_hectare_bordas.txt'
	try:
		viz_file=open(arquivoviz,"rb")
		viz_pck=viz_file.read()
		viz_file.close()
		lista_borda=pck.loads(viz_pck)
	except FileNotFoundError:
		print("Não achou a lista de patches na borda. Iniciando processo de criação.")
		lista_borda=criando_lista_patches_cada_borda()
		a=open(arquivoviz,"ab")		 
		dados_em_pck=pck.dumps(lista_borda)   
		a.write(dados_em_pck)
		a.close()
	#aqui vai checar lista total
	arquivoviz=f'50x50_patches_por_hectare.txt'
	try:
		viz_file=open(arquivoviz,"rb")
		viz_pck=viz_file.read()
		viz_file.close()
		lista_patches=pck.loads(viz_pck)
	except FileNotFoundError:
		print("Não achou a lista de patches do hectare. Iniciando processo de criação.")
		lista_patches=criando_lista_patches_cada_hectare()
		a=open(arquivoviz,"ab")		 
		dados_em_pck=pck.dumps(lista_patches)   
		a.write(dados_em_pck)
		a.close()
	return viz,lista_patches,lista_borda
def criando_lista_viz(lin,col,distancia_maxima):
	lista_raios=[]
	viz=[] #criando uma matriz vazia
	###bloco de criar a matriz de vizinhos vazia para poder adicionar os vizinhos
	#e como ja vai percorrer uma vez a matriz inteira, eu aproveito e ja calculo todos os vizinhos possiveis e faco a lista
	#de vizinhos possiveis.
	for l in range(lin): #percorrendo os loops de linhas
		viz.append([]) #colocando as linhas na matriz
		for c in range(col): #percorrendo os loops de colunas
			viz[l].append([]) #colocando as colunas na matriz
			if l==0 and c==0: #aqui checo soh pra ver se eh o primeiro patch de todos
				continue #se for ele nao faz a conta, pois ele eh raio zero
			else:	#qnd n for o patch 0,0 ele vai entrar aqui
				distancia_tempo=((0-l)**2+(0-c)**2)**(1/2)
				if distancia_tempo<=distancia_maxima:
					lista_raios.append(distancia_tempo)

					 #aqui eu faco a conta da distancia e adiciona na lista
	lista_raios=sorted(set(lista_raios)) #qnd o loop termina, eu retiro os valores repetidos e ordeno a lista
	#pois agora ela vai ter todas as distancias presentes nesse grid, e ordenadas por raio.

	###bloco de adicionar listas dos raios possiveis dentro da matriz
	for l in range(lin):
		for c in range(col):
			for i in range(len(lista_raios)):
				viz[l][c].append([]) #colocando todos os raios possiveis da matriz

	###bloco de achar os vizinhos e por no raio correspondente
	for l in range(lin): #loop da linha para o patch focal
		for c in range(col): #loop da coluna para o patch focal
			for ll in range(-lin,lin): #loop de linha do vizinho do patch focal
				for cc in range(-col,col): #loop de coluna do vizinho do patch focal
					if l==ll and c==cc: #checando se o patch focal eh o mesmo q o vizinho
						continue #se for, acontece nada
					else: 
						dist=((l-ll)**2+(c-cc)**2)**(1/2) #calculando a distancia do patch focal(l,c) pro vizinho(ll,cc)
						for i,j in enumerate(lista_raios): #loop para rodar nos raios existentes, pegando a casa q ele esta(i) e o valor do raio (j)
							if dist==j: #checando se o valor da distancia calculado eh igual a algum raio e descobrindo em qual lista de vizinho por
								if ll>=lin:
									tempolin=(lin-1)+((lin-1)-ll)
								elif ll<0:
									tempolin=ll*(-1)
								else:
									tempolin=ll
								if cc>=col:
									tempocol=(col-1)+((col-1)-cc)
								elif cc<0:
									tempocol=cc*(-1)
								else:
									tempocol=cc
								viz[l][c][i].append([tempolin,tempocol])
								break										
	return viz,lista_raios
def criando_popmigrante_densidade_patch(l,c,t,H0,h0,f0,fracao_indv_migrante,g_all,pop_migrante,emigracao_patchs):
	#p
	if g_all[l][c][2][t]==0:
		pop_migrante[l][c][0].append(int(g_all[l][c][0][t]))
		emigracao_patchs[l][c][0].append(int(g_all[l][c][0][t])) 
	else:
		pop_migrante[l][c][0].append(int(fracao_indv_migrante[0]*(H0[0]/(H0[0]+g_all[l][c][2][t]))*((g_all[l][c][0][t]**2)/(g_all[l][c][0][t]+f0[0]))))
		emigracao_patchs[l][c][0].append(int(fracao_indv_migrante[0]*(H0[0]/(H0[0]+g_all[l][c][2][t]))*((g_all[l][c][0][t]**2)/(g_all[l][c][0][t]+f0[0]))))
	#q
	if g_all[l][c][2][t]==0:#and g_all[l][c][3][t]==0:
		pop_migrante[l][c][1].append(int(g_all[l][c][1][t]))
		emigracao_patchs[l][c][1].append(int(g_all[l][c][1][t]))
	else:

		pop_migrante[l][c][1].append(int(fracao_indv_migrante[1]*(H0[1]/(H0[1]+(g_all[l][c][2][t]+g_all[l][c][3][t])))*((g_all[l][c][1][t]**2)/(g_all[l][c][1][t]+f0[1]))))
		emigracao_patchs[l][c][1].append(int(fracao_indv_migrante[1]*(H0[1]/(H0[1]+(g_all[l][c][2][t]+g_all[l][c][3][t])))*((g_all[l][c][1][t]**2)/(g_all[l][c][1][t]+f0[1]))))
	#n
	pop_migrante[l][c][2].append(int((fracao_indv_migrante[2]*(g_all[l][c][2][t]**2))/(g_all[l][c][2][t]+h0[0])))
	emigracao_patchs[l][c][2].append(int((fracao_indv_migrante[2]*(g_all[l][c][2][t]**2))/(g_all[l][c][2][t]+h0[0])))
	#h 
	pop_migrante[l][c][3].append(int((fracao_indv_migrante[3]*(g_all[l][c][3][t]**2))/(g_all[l][c][3][t]+h0[1])))
	emigracao_patchs[l][c][3].append(int((fracao_indv_migrante[3]*(g_all[l][c][3][t]**2))/(g_all[l][c][3][t]+h0[1])))

	return pop_migrante,emigracao_patchs
def sorteio_vizinhos(g_all,pop_migrante,taxa_disp,viz,l,c,t,visitacao_patchs,imigracao_patchs,lista_distancia):
	#criando a pop migrante temporaria para poder diminuir desta variavel ateh chegar a zero
	pop_migrante_tempo=[]
	for i,p in enumerate(taxa_disp):
		pop_migrante_tempo.append(pop_migrante[l][c][i][t])
	for i,p in enumerate(pop_migrante_tempo): #i eh a casa. p eh o elemento #loop para rodar a qnt certa de pop e ja usar as taxa dela de dispersao
		x=0
		while pop_migrante_tempo[i]>0 and x<len(lista_distancia):
			for j in viz[l][c][x]:
				
				sorteado=random.choice(viz[l][c][x])
				indv=(taxa_disp[i]/lista_distancia[x])*pop_migrante[l][c][i][t]
				if indv<1: #conferindo se vai migrar um valor menor do que 1, se for eu arredondo para 1
					indv=1
				if indv>pop_migrante_tempo[i]:
					indv=pop_migrante_tempo[i]
				pop_migrante_tempo[i]-=indv
				indv=int(indv)
				g_all[sorteado[0]][sorteado[1]][i][t]+=indv
				visitacao_patchs[sorteado[0]][sorteado[1]][i][t]+=1
				imigracao_patchs[sorteado[0]][sorteado[1]][i][t]+=indv
				if pop_migrante_tempo[i]<=0:
					break
			x+=1
	return g_all,visitacao_patchs,imigracao_patchs
def dinamica_tempo0(g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante,c,l,t,patch_iniciais):
	achou_patch=0
	for i,j in enumerate(patch_iniciais): #aqui eu fiz um for para passar pela lista de patchs iniciais para
		if j[0]==l and j[1]==c: #eu poder comparar com qual linha/coluna está. se alguma linha/coluna for igual ao da lista
			#de patchs iniciais, ele adiciona a pop inicial. o 'i' serve para saber qual casa o patch esta na lista patch_iniciais
			#e ai a mesma casa é usada para achar a pop inicial daquele patch.
			
			#pop 0 - P
			achou_patch=1
			g_all[l][c][0].append(pop_iniciais[i][0])
			if g_all[l][c][0][t]>0.01:
				ocupacao_total_tempo[0]+=1
				media_regional_tempo[0]+=g_all[l][c][0][t]


			#pop 1 - Q
			#testando
			g_all[l][c][1].append(pop_iniciais[i][1])
			if g_all[l][c][1][t]>0.01:
				ocupacao_total_tempo[1]+=1
				media_regional_tempo[1]+=g_all[l][c][1][t]


			#pop 2 - N
			g_all[l][c][2].append(pop_iniciais[i][2])
			if g_all[l][c][2][t]>0.01:
				ocupacao_total_tempo[2]+=1
				media_regional_tempo[2]+=g_all[l][c][2][t]
			
			
			#pop 3 - H
			g_all[l][c][3].append(pop_iniciais[i][3])
			if g_all[l][c][3][t]>0.01:
				ocupacao_total_tempo[3]+=1
				media_regional_tempo[3]+=g_all[l][c][3][t]


			#colocando zero na migracao
			pop_migrante[l][c][0].append(0)
			pop_migrante[l][c][1].append(0)
			pop_migrante[l][c][2].append(0)
			pop_migrante[l][c][3].append(0)

	if achou_patch==0:
		#colocando zero nas pop iniciais
		g_all[l][c][0].append(0)
		g_all[l][c][1].append(0)
		g_all[l][c][2].append(0)
		g_all[l][c][3].append(0)

		#colocando zero na migracao
		pop_migrante[l][c][0].append(0)
		pop_migrante[l][c][1].append(0)
		pop_migrante[l][c][2].append(0)
		pop_migrante[l][c][3].append(0)
	return g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante
def dinamica_iterativa(g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante,c,l,t,a1,a2,a3,lambda1,lambda2,TS,TH,H0,h0,f0,ocupacao_total,Pnumerico,Q1numerico,Q2numerico):
	
	rfP=(math.exp((-a1*TS[0]*g_all[l][c][0][t-1])/(1+a1*TH[0]*g_all[l][c][2][t-1])))#resposta funcional especialista(P) (P para o N)
	rfQ1=(math.exp((-a2*TS[1]*g_all[l][c][1][t-1])/(1+a2*TH[1]*g_all[l][c][2][t-1])))#resposta funcional generalsita(Q) (Q para o N)
	rfQ2=(math.exp((-a3*TS[2]*g_all[l][c][1][t-1])/(1+a3*TH[2]*g_all[l][c][3][t-1])))#resposta funcional generalsita(Q) (Q para o H) #to-do: conferir a resposta funcional do generalista com o host alternativo (q->h)
	
	#vou checar se a resposta funcional é zero ou menor, pq se for a equação dos bixos daria zero pq multiplica ela
	if rfP<=0:
		rfP=1 #coloco igual a 1 pq ai a resposta funcional nao vai alterar em nada as equações.
	if rfQ1<=0:
		rfQ1=1
	if rfQ2<=0:
		rfQ2=1

	#bloquinho dinamica populacao 0 (P)
	if g_all[l][c][0][t-1]==0:
		g_all[l][c][0].append(0)
	else:		
		g_all[l][c][0].append(g_all[l][c][2][t-1]*(1-rfP)*Pnumerico)
		if g_all[l][c][0][t]<=0.01:
			g_all[l][c][0][t]=0

	#bloquinho dinamica populacao 1 (Q):
	#teste
	if g_all[l][c][1][t-1]==0:
		g_all[l][c][1].append(0)
	else:
		g_all[l][c][1].append(g_all[l][c][2][t-1]*rfP*(1-rfQ1)*Q1numerico) #to-do: essa equação só leva em conta 1 hospedeiro.
		if g_all[l][c][1][t]<=0.01:
			g_all[l][c][1][t]=0

	#bloquinho dinamica populacao 2 (N)
	if g_all[l][c][2][t-1]==0:
		g_all[l][c][2].append(0)
	else:
		g_all[l][c][2].append(g_all[l][c][2][t-1]*lambda1*rfP*rfQ1)
		if g_all[l][c][2][t]<=0.01:
			g_all[l][c][2][t]=0
		if g_all[l][c][2][t]>1000 and ocupacao_total[0][t-1]==0 and ocupacao_total[1][t-1]==0: 
			g_all[l][c][2][t]=1000


	#bloquinho dinamica populacao 3 (h):
	#teste:
	if g_all[l][c][3][t-1]==0:
		g_all[l][c][3].append(0)
	else:
		g_all[l][c][3].append(g_all[l][c][3][t-1]*lambda2*rfQ2) #to-do: conferir a equação do host alternativo (h)
		if g_all[l][c][3][t]<=0.01:
			g_all[l][c][3][t]=0
		if g_all[l][c][3][t]>1000 and ocupacao_total[1][t-1]==0:
			g_all[l][c][3][t]=1000
	return g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante
def conferir_recurso_parasitoides(g_all,l,c,t):
	if g_all[l][c][0][t]>0 and g_all[l][c][2][t]==0:
		g_all[l][c][0][t]=0
	if g_all[l][c][1][t]>0 and g_all[l][c][2][t]==0 and g_all[l][c][3][t]==0:
		g_all[l][c][1][t]=0
	return g_all
def salva_arquivos(zzzz,g_all,pop_migrante,ocupacao_total,media_regional,media_regional_migracao,visitacao_patchs,imigracao_patchs,emigracao_patchs,a1,a2,a3,lambda1,lambda2,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,qnd_onde_liberou):
	

	nome=f'{zzzz}_qnd_onde_liberou.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(qnd_onde_liberou)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_g_all.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(g_all)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_pop_migrante.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(pop_migrante)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_ocupacao_total.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(ocupacao_total)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_media_regional.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(media_regional)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_media_regional_migracao.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(media_regional_migracao)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_visitacao_patchs.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(visitacao_patchs)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_imigracao_patchs.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(imigracao_patchs)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_emigracao_patchs.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(emigracao_patchs)   
	a.write(dados_em_pck)
	a.close()

	nome=f'read_me.txt'
	a=open(nome,"w")
	conteudo_readme=f'A ordem do nome é tipo.do.arquivo.txt\n \n \n Os arquivos possuem listas:\n g_all[linha][coluna][populacao][tempo]: é a populacao de cada patch por tempo\npop_migrante[linha][coluna][populacao][tempo]: é o vetor de migracao\nocupacao_total[populacao][tempo]: é a ocupacao de patchs por tempo\nmedia_regional[populacao][tempo]: é a media regional de populacao \nmedia_regional_migracao[populacao][tempo]: é a media regional de migracao\nvisitacao_patchs[linha][coluna][populacao][tempo]:é quais patchs receberam imigrantes\nimigracao_patchs[linha][coluna][populacao][tempo]: é o quanto cada patch recebeu de migrantes, ou seja, a imigracao\nemigracao_patchs[linha][coluna][pop][tempo]:é o qnt que cada pop de cada patch contribuiu para a migraçao, ou seja, a emigraçao \n \nOs valores desse conjunto de simulação são:\na1={a1}\na2={a2}\na3={a3}\nlambda1={lambda1}\nlambda2={lambda2}\ntaxa_disp={taxa_disp}\nfracao_indv_migrante={fracao_indv_migrante}\ntempo_final={tempo_final}\nlin={lin}\ncol={col}'
	a.write(conteudo_readme)
	a.close()
	return
def criador_vetor(lin,col,tempo_final):
	#Primeiro eu crio uma lista vazia, e dps percorro o numero de linhas e vou adicionando uma lista vazia opr iteracao.
	#Em seguida, dentro de cada linha, eu percorro o numero de colunas e vou adicionando lista vazia por coluna, por linha
	#e em seguida percorro o numero de pop dentro de cada coluna e adiciona listas vazias iguais o numero de pop
	#dentro de cada coluna dentro de cada linha.
	g_all=[] #g_all[lin][col][pop][tempo]
	pop_migrante=[] #pop_migrante[lin][col][pop][tempo]
	visitacao_patchs=[]
	imigracao_patchs=[]
	emigracao_patchs=[]
	for l in range(lin):
		g_all.append([])
		pop_migrante.append([])
		visitacao_patchs.append([])
		imigracao_patchs.append([])
		emigracao_patchs.append([])
		for c in range(col):
			g_all[l].append([[],[],[],[]])
			pop_migrante[l].append([[],[],[],[]])
			visitacao_patchs[l].append([[],[],[],[]])
			imigracao_patchs[l].append([[],[],[],[]])
			emigracao_patchs[l].append([[0],[0],[0],[0]])
			for t in range(tempo_final):
				visitacao_patchs[l][c][0].append(0)
				visitacao_patchs[l][c][1].append(0)
				visitacao_patchs[l][c][2].append(0)
				visitacao_patchs[l][c][3].append(0)

				imigracao_patchs[l][c][0].append(0)
				imigracao_patchs[l][c][1].append(0)
				imigracao_patchs[l][c][2].append(0)
				imigracao_patchs[l][c][3].append(0)
	return g_all,pop_migrante,visitacao_patchs,imigracao_patchs,emigracao_patchs
def sorteio_pop_inicial_host(lin,col,porcentagem_popinihost,inicial_h,inicial_n,solt_valor,loop,lista_patches): #retorna em primeiro os patches iniciais e em segundo as populações iniciais APENAS DO HOST
	pop_iniciais=[]
	grid=lin*col
	numero_patch_lib=grid*porcentagem_popinihost[0]
	numero_patch_lib/=10
	if loop==21 or loop==22 or loop==25 or loop==26 or loop==27 or loop==28:
		porcentagem_n=int(numero_patch_lib-(solt_valor*2)) #n
		porcentagem_h=int(numero_patch_lib-(solt_valor*2)) #h #aqui ta indo com os valores do outro host
	else:
		porcentagem_n=int(numero_patch_lib-solt_valor) #n
		porcentagem_h=int(numero_patch_lib-solt_valor) #h #aqui ta indo com os valores do outro host
	patch_ini=[]

	for posicao_patch,lista_dos_patches in enumerate(lista_patches):
		for i in range(porcentagem_n):
			patch_tempo=random.choice(lista_dos_patches)
			while patch_tempo in patch_ini:
				patch_tempo=random.choice(lista_dos_patches)
			patch_ini.append(patch_tempo)
			pop_iniciais.append([0,0,inicial_n,0])
	return patch_ini,pop_iniciais
def alocar_parasitoide(borda_meio,solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,hec10_bordas,inicial_n,inicial_h,para_lista_inicial,loop):
	if borda_meio==0:
		patches_sorteados=[[],[]]
		tempo=BORDA_alocando_parasitoide(solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,lista_borda,inicial_n,inicial_h,para_lista_inicial,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		return patch_iniciais,pop_iniciais,patches_sorteados

	else:
		tempo=MEIO_alocando_parasitoide(solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,inicial_n,inicial_h,para_lista_inicial,loop,lista_patches)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		patches_sorteados=tempo[2]
		return patch_iniciais,pop_iniciais,patches_sorteados
def MEIO_alocando_parasitoide(solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,inicial_n,inicial_h,para_lista_inicial,loop,lista_patches):
	if loop==19 or loop==20 or loop==23 or loop==24:
		patches_sorteados=[[],[]]
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],para_lista_inicial[1][solt_posi],inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
	elif loop==21 or loop==25:
		patches_sorteados=[[],[]]
		#pondo cotesia:
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo tetrastichus
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([0,para_lista_inicial[1][solt_posi],inicial_n,0])
				patches_sorteados[1].append(patch_tempo)
	elif loop==22 or loop==26 or loop==27:
		patches_sorteados=[[],[]]
		#pondo cotesia:
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo segunda cotesia
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
	elif loop==28:
		patches_sorteados=[[],[]]
		#pondo tetrastichus:
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([0,para_lista_inicial[1][solt_posi],inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo segunda tetrastichus
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([0,para_lista_inicial[1][solt_posi],inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
	
	return patch_iniciais,pop_iniciais,patches_sorteados
def BORDA_alocando_parasitoide(solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,lista_borda,inicial_n,inicial_h,para_lista_inicial,loop):
	
	if loop==19 or loop==20 or loop==23 or loop==24:
		patches_sorteados=[[],[]]
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],para_lista_inicial[1][solt_posi],inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
	elif loop==21 or loop==25:
		patches_sorteados=[[],[]]
		#pondo cotesia:
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo tetrastichus
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([0,para_lista_inicial[1][solt_posi],inicial_n,0])
				patches_sorteados[1].append(patch_tempo)
	elif loop==22 or loop==26 or loop==27:
		patches_sorteados=[[],[]]
		#pondo cotesia:
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo segunda cotesia
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
	elif loop==28:
		patches_sorteados=[[],[]]
		#pondo tetrastichus:
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([0,para_lista_inicial[1][solt_posi],inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo segunda tetrastichus
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([0,para_lista_inicial[1][solt_posi],inicial_n,0])
				patches_sorteados[0].append(patch_tempo)
	
	return patch_iniciais,pop_iniciais,patches_sorteados
def checando_necessidade_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,lista_patches,lista_borda,qnd_onde_liberou):
	#lista_patches,lista_borda
	teve_liberacao=False
	pop_n_inicial=0
	pop_n_t=0
	for posicao_hectare,hectare in enumerate(lista_patches):
		pop_n_inicial=0
		pop_n_t=0
		for patch in hectare:
			pop_n_t+=g_all[patch[0]][patch[1]][2][t]
		if pop_n_t>=925: #se mudar a pop inicial tem q mudar isso.
			teve_liberacao=True
			qnd_onde_liberou.append(f't{t}h{posicao_hectare}')
			if borda_meio==0: #borda
				patches_para_liberar=lista_borda[posicao_hectare]
				g_all=nova_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,patches_para_liberar)
			else: #meio
				patches_para_liberar=lista_patches[posicao_hectare]
				g_all=nova_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,patches_para_liberar)
	return g_all,qnd_onde_liberou,teve_liberacao
def nova_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,patches_para_liberar):
	#aqui eu preciso garantir q qnd é separado eu GARANTA q nao acontece de ficar igual
	if loop==21 or loop==19:
		patches_q_um_para_foi=[]
		#pondo cotesia
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][0][t]+=para_lista_inicial[0][solt_posi]
			patches_q_um_para_foi.append(patch_tempo)
		#pondo tetrastichus:
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][1][t]+=para_lista_inicial[1][solt_posi]
	elif loop==20 : 
		patches_q_um_para_foi=[]
		#pondo ambos
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][0][t]+=para_lista_inicial[0][solt_posi]
			g_all[patch_tempo[0]][patch_tempo[1]][1][t]+=para_lista_inicial[1][solt_posi]
			patches_q_um_para_foi.append(patch_tempo)
	elif loop==22:
		patches_q_um_para_foi=[]
		#pondo primeira cotesia
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][0][t]+=para_lista_inicial[0][solt_posi]
			patches_q_um_para_foi.append(patch_tempo)
		#pondo segunda cotesia:
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][0][t]+=para_lista_inicial[0][solt_posi]

	return g_all
def loops_cenarios(loop,viz,lista_patches,lista_borda):
	if loop>5: #solturas qnd tem 10 hectares 
		#10 hectares:
		#1 no meio
		hec10_1_meio_separados1=[[4,7],[4,32],[14,7],[14,32],[24,7],[24,32],[34,7],[34,32],[44,7],[44,32]]#[[4,7],[4,32],[14,7],[14,32],[24,7],[24,32],[34,7],[34,32],[44,7],[44,32]]
		hec10_1_meio_separados2=[[4,17],[4,42],[14,17],[14,42],[24,17],[24,42],[34,17],[34,42],[44,17],[44,42]]#[[4,17],[4,42],[14,17],[14,42],[24,17],[24,42],[34,17],[34,42],[44,17],[44,42]]

		#1 na borda
		hec10_1_borda_separados1=[[0,7],[0,32],[13,0],[13,49],[23,0],[33,0],[23,49],[33,49],[49,7],[49,32]]
		hec10_1_borda_separados2=[[0,17],[0,42],[17,0],[17,49],[27,0],[37,0],[27,49],[37,49],[49,17],[49,42]]

		#4 na borda
		hec10_4_borda_separados1=[[2,0],[0,9],[0,19],[7,23],[7,24],[0,33],[0,43],[2,49],[42,0],[49,9],[49,19],[47,23],[47,24],[49,33],[49,43],[42,49],[13,0],[19,0],[10,23],[16,23],[10,24],[16,24],[13,49],[19,49],[23,0],[29,0],[20,23],[26,23],[20,24],[26,24],[23,49],[29,49],[33,0],[39,0],[30,23],[36,23],[30,24],[36,24],[33,49],[39,49]]
		hec10_4_borda_separados2=[[7,0],[0,4],[0,14],[2,23],[2,24],[0,28],[0,38],[7,49],[47,0],[49,4],[49,14],[42,23],[42,24],[49,28],[49,38],[47,49],[10,0],[16,0],[13,23],[19,23],[13,24],[19,24],[10,49],[16,49],[20,0],[26,0],[23,23],[29,23],[23,24],[29,24],[20,49],[26,49],[30,0],[36,0],[33,23],[39,23],[33,24],[39,24],[30,49],[36,49]]

		#4 no meio:
		hec10_4_meio_separados1=[[2,9],[2,19],[7,9],[7,19],[2,34],[2,44],[7,34],[7,44],[12,9],[12,19],[17,9],[17,19],[12,34],[12,44],[17,34],[17,44],[22,9],[22,19],[27,9],[27,19],[22,34],[22,44],[27,34],[27,44],[32,9],[32,19],[37,9],[37,19],[32,34],[32,44],[37,34],[37,44],[42,9],[42,19],[47,9],[47,19],[42,34],[42,44],[47,34],[47,44]]
		hec10_4_meio_separados2=[[2,4],[2,14],[7,4],[7,14],[2,29],[2,39],[7,29],[7,39],[12,4],[12,14],[17,4],[17,14],[12,29],[12,39],[17,29],[17,39],[22,4],[22,14],[27,4],[27,14],[22,29],[22,39],[27,29],[27,39],[32,4],[32,14],[37,4],[37,14],[32,29],[32,39],[37,29],[37,39],[42,4],[42,14],[47,4],[47,14],[42,29],[42,39],[47,29],[47,39]]
		
		#6 no meio
		#hec10_4_meio_separados1=[[2,2],[7,2],[2,10],[7,10],[2,18],[7,18],[12,2],[17,2],[12,10],[17,10],[12,18],[17,18],[22,2],[27,2],[22,10],[7,10],[22,18],[27,18],[32,2],[37,2],[32,10],[37,10],[32,18],[37,18],[42,2],[47,2],[42,10],[47,10],[42,18],[47,18],[2,27],[7,27],[2,35],[7,35],[2,43],[7,43],[12,27],[17,27],[12,35],[17,35],[12,43],[17,43],[22,27],[27,27],[22,35],[7,35],[22,43],[27,43],[32,27],[37,27],[32,35],[37,35],[32,43],[37,43],[42,27],[47,27],[42,35],[47,35],[42,43],[47,43]]
		#hec10_4_meio_separados2=[[2,6],[7,6],[2,14],[7,14],[2,22],[7,22],[12,6],[17,6],[12,14],[17,14],[12,22],[17,22],[22,6],[27,6],[22,14],[27,14],[22,22],[27,22],[32,6],[37,6],[32,14],[37,14],[32,22],[37,22],[42,6],[47,6],[42,14],[47,14],[42,22],[47,22],[2,31],[7,31],[2,39],[7,39],[2,47],[7,47],[12,31],[17,31],[12,39],[17,39],[12,47],[17,47],[22,31],[27,31],[22,39],[27,39],[22,47],[27,47],[32,31],[37,31],[32,39],[37,39],[32,47],[37,47],[42,31],[47,31],[42,39],[47,39],[42,47],[47,47]]
		#4 na borda
		#hec10_4_borda_separados1=[[1,0],[6,0],[0,2],[0,12],[1,23],[6,23],[1,24],[6,24],[0,27],[0,37],[1,49],[6,49],[48,0],[43,0],[49,2],[49,12],[48,23],[43,23],[48,24],[43,24],[49,27],[49,37],[48,49],[43,49],[11,0],[14,0],[17,0],[12,23],[15,23],[18,23],[12,24],[15,24],[18,24],[11,49],[14,49],[17,49],[21,0],[24,0],[27,0],[22,23],[25,23],[28,23],[22,24],[25,24],[28,24],[21,49],[24,49],[27,49],[31,0],[34,0],[37,0],[32,23],[35,23],[38,23],[32,24],[35,24],[38,24],[31,49],[34,49],[37,49]]
		#hec10_4_borda_separados2=[[3,0],[8,0],[0,6],[0,18],[3,23],[8,23],[3,24],[8,24],[0,31],[0,43],[3,49],[8,49],[46,0],[41,0],[49,6],[49,18],[46,23],[41,23],[46,24],[41,24],[49,31],[49,43],[46,49],[41,49],[12,0],[15,0],[18,0],[11,23],[14,23],[17,23],[11,24],[14,24],[17,24],[12,49],[15,49],[18,49],[22,0],[25,0],[28,0],[21,23],[24,23],[27,23],[21,24],[24,24],[27,24],[22,49],[25,49],[28,49],[32,0],[35,0],[38,0],[31,23],[34,23],[37,23],[31,24],[34,24],[37,24],[32,49],[35,49],[38,49]]

		#8 no meio
		hec10_8_meio_separados1=[[2,1],[7,1],[2,7],[7,7],[2,13],[7,13],[2,19],[7,19],[12,1],[17,1],[12,7],[17,7],[12,13],[17,13],[12,19],[17,19],[22,1],[27,1],[22,7],[27,7],[22,13],[27,13],[22,19],[27,19],[32,1],[37,1],[32,7],[37,7],[32,13],[37,13],[32,19],[37,19],[42,1],[47,1],[42,7],[47,7],[42,13],[47,13],[42,19],[47,19],[2,26],[7,26],[2,32],[7,32],[2,38],[7,38],[2,44],[7,44],[12,26],[17,26],[12,32],[17,32],[12,38],[17,38],[12,44],[17,44],[22,26],[27,26],[22,32],[27,32],[22,38],[27,38],[22,44],[27,44],[32,26],[37,26],[32,32],[37,32],[32,38],[37,38],[32,44],[37,44],[42,26],[47,26],[42,32],[47,32],[42,38],[47,38],[42,44],[47,44]]
		hec10_8_meio_separados2=[[2,4],[7,4],[2,10],[7,10],[2,16],[7,16],[2,22],[7,22],[12,4],[17,4],[12,10],[17,10],[12,16],[17,16],[12,22],[17,22],[22,4],[27,4],[22,10],[27,10],[22,16],[27,16],[22,22],[27,22],[32,4],[37,4],[32,10],[37,10],[32,16],[37,16],[32,22],[37,22],[42,4],[47,4],[42,10],[47,10],[42,16],[47,16],[42,22],[47,22],[2,29],[7,29],[2,35],[7,35],[2,41],[7,41],[2,47],[7,47],[12,29],[17,29],[12,35],[17,35],[12,41],[17,41],[12,47],[17,47],[22,29],[27,29],[22,35],[27,35],[22,41],[27,41],[22,47],[27,47],[32,29],[37,29],[32,35],[37,35],[32,41],[37,41],[32,47],[37,47],[42,29],[47,29],[42,35],[47,35],[42,41],[47,41],[42,47],[47,47]]
		#8 na borda
		hec10_8_borda_separados1=[[1,0],[6,0],[0,2],[0,12],[1,23],[6,23],[0,9],[0,21],[1,24],[6,24],[0,27],[0,37],[1,49],[6,49],[0,34],[0,47],[48,0],[43,0],[49,2],[49,12],[48,23],[43,23],[49,9],[49,21],[48,24],[43,24],[49,27],[49,37],[48,49],[43,49],[49,34],[49,47],[11,0],[14,0],[17,0],[12,23],[15,23],[18,23],[16,0],[13,23],[12,24],[15,24],[18,24],[11,49],[14,49],[17,49],[13,24],[16,49],[21,0],[24,0],[27,0],[22,23],[25,23],[28,23],[26,0],[23,23],[22,24],[25,24],[28,24],[21,49],[24,49],[27,49],[23,24],[26,49],[31,0],[34,0],[37,0],[32,23],[35,23],[38,23],[36,0],[33,23],[32,24],[35,24],[38,24],[31,49],[34,49],[37,49],[33,24],[36,49]]
		hec10_8_borda_separados2=[[3,0],[8,0],[0,6],[0,18],[3,23],[8,23],[0,4],[0,15],[3,24],[8,24],[0,31],[0,43],[3,49],[8,49],[0,29],[0,40],[46,0],[41,0],[49,6],[49,18],[46,23],[41,23],[49,4],[49,15],[46,24],[41,24],[49,31],[49,43],[46,49],[41,49],[49,29],[49,40],[12,0],[15,0],[18,0],[11,23],[14,23],[17,23],[13,0],[16,23],[11,24],[14,24],[17,24],[12,49],[15,49],[18,49],[16,24],[13,49],[22,0],[25,0],[28,0],[21,23],[24,23],[27,23],[23,0],[26,23],[21,24],[24,24],[27,24],[22,49],[25,49],[28,49],[26,24],[23,49],[32,0],[35,0],[38,0],[31,23],[34,23],[37,23],[33,0],[36,23],[31,24],[34,24],[37,24],[32,49],[35,49],[38,49],[36,24],[33,49]]
		
		hec10_8_meio=[hec10_8_meio_separados1,hec10_8_meio_separados2]
		hec10_8_borda=[hec10_8_borda_separados1,hec10_8_borda_separados2]
		hec10_4_borda=[hec10_4_borda_separados1,hec10_4_borda_separados2]
		hec10_4_meio=[hec10_4_meio_separados1,hec10_4_meio_separados2]
		hec10_1_borda=[hec10_1_borda_separados1,hec10_1_borda_separados2]
		hec10_1_meio=[hec10_1_meio_separados1,hec10_1_meio_separados2]
		
		

		hec10_1=[hec10_1_borda,hec10_1_meio]#sempre primeiro a borda e dps o meio
		hec10_4=[hec10_4_borda,hec10_4_meio]#sempre primeiro a borda e dps o meio
		hec10_8=[hec10_8_borda,hec10_8_meio]#sempre primeiro a borda e dps o meio
		hec10=[hec10_1,hec10_4,hec10_8]
		hec10_bordas=[hec10_1_borda,hec10_4_borda,hec10_8_borda]
	elif loop>19 and loop>4: #solturas qnd tem 5 hectares
		solturas5pontos=[[[25,5],[25,15],[25,25],[25,35],[25,45]],[[25,25],[0,25],[49,25],[25,0],[25,49]]]
		solturas1ponto=[[25,25],[0,25]]
	##############################
	####CENARIOS COM LIBERAÇÃO####
	##############################
	#liberando 6k cotesia e 7k tetrastichus por hectare
	if loop==23: #só 6k cotesia
		Pnumerico=52.5
		Q1numerico=30
		Q2numerico=0
		replicas=0 #replicas
		TS=[4,4,0] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,0]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=0
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=0
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.25,0,0.05,0]#[0.125,0.125,0.05,0]  #PQNH
		tempo_final=31
		inicial_n=41
		inicial_h=0
		#
		solturas=[8]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[750],[0]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[[0.9,0.9,0.3,0]]#taxa de qnts parasitoides saem do patch baixa e alta. 
		porcentagens=[0.15] #ocupacao inicial do host
		for solt_posi,solt_valor in enumerate(solturas):
			for pop_host_ini_porcentagem in porcentagens:
				porcentagem_popinihost=[pop_host_ini_porcentagem,0]
				for fracao_indv_migrante in indv_migrante:
					for borda_meio,hec10_BordaMeio in enumerate(hec10[solt_posi]): #pega borda e meio
						
						tempo=sorteio_pop_inicial_host(lin,col,porcentagem_popinihost,inicial_h,inicial_n,solt_valor,loop,lista_patches)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]

						
						tempo=alocar_parasitoide(borda_meio,solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,hec10_bordas,inicial_n,inicial_h,para_lista_inicial,loop)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]	
						patches_sorteados=tempo[2]				
						if borda_meio==0:
							bOUm="borda"
						else:
							bOUm="meio"
						wd_novo=f'{wd_original}\\{bOUm}'
						os.makedirs(wd_novo)
						os.chdir(wd_novo)
						zzzz=0
						replicas=0
						total_replicas=0
						while replicas<10 or total_replicas<50:
							modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,patches_sorteados,bOUm,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,borda_meio,loop,lib_extras,Pnumerico,Q1numerico,Q2numerico)
							ocupacao_file=open(f'{zzzz}_ocupacao_total.txt',"rb")
							zzzz+=1
							ocupacao_pck=ocupacao_file.read()
							ocupacao_file.close()
							ocupacao=pck.loads(ocupacao_pck)
							total_replicas+=1
							if total_replicas>49:
								break
							if ocupacao[0][30]>0:
								replicas+=1		
	elif loop==24: #ambos parasitoides juntos 6k cotesia e 7k tetrastichus
		Pnumerico=52.5
		Q1numerico=30
		Q2numerico=0
		replicas=0 #replicas
		TS=[4,4,0] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,0]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=0
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=0
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.25,0.25,0.05,0]#[0.125,0.125,0.05,0]  #PQNH
		tempo_final=31
		inicial_n=41
		inicial_h=0
		#
		solturas=[8]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[750],[875]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[[0.9,0.9,0.3,0]]#taxa de qnts parasitoides saem do patch baixa e alta. 
		porcentagens=[0.15] #ocupacao inicial do host
		for solt_posi,solt_valor in enumerate(solturas):
			for pop_host_ini_porcentagem in porcentagens:
				porcentagem_popinihost=[pop_host_ini_porcentagem,0]
				for fracao_indv_migrante in indv_migrante:
					for borda_meio,hec10_BordaMeio in enumerate(hec10[solt_posi]): #pega borda e meio
						
						tempo=sorteio_pop_inicial_host(lin,col,porcentagem_popinihost,inicial_h,inicial_n,solt_valor,loop,lista_patches)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]
						
						tempo=alocar_parasitoide(borda_meio,solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,hec10_bordas,inicial_n,inicial_h,para_lista_inicial,loop)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]	
						patches_sorteados=tempo[2]				
						if borda_meio==0:
							bOUm="borda"
						else:
							bOUm="meio"
						wd_novo=f'{wd_original}\\{bOUm}'
						os.makedirs(wd_novo)
						os.chdir(wd_novo)
						zzzz=0
						replicas=0
						total_replicas=0
						while replicas<10 or total_replicas<50:
							modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,patches_sorteados,bOUm,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,borda_meio,loop,lib_extras,Pnumerico,Q1numerico,Q2numerico)

							ocupacao_file=open(f'{zzzz}_ocupacao_total.txt',"rb")
							zzzz+=1
							ocupacao_pck=ocupacao_file.read()
							ocupacao_file.close()
							ocupacao=pck.loads(ocupacao_pck)
							total_replicas+=1
							if total_replicas>49:
								break
							if ocupacao[0][30]>0:
								replicas+=1
	elif loop==25: #ambos parasitoides separados 6k cotesia e 7k tetrastichus
		Pnumerico=52.5
		Q1numerico=30
		Q2numerico=0
		replicas=0 #replicas
		TS=[4,4,0] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,0]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=0
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=0
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.25,0.25,0.05,0]#[0.125,0.125,0.05,0]  #PQNH
		tempo_final=31
		inicial_n=41
		inicial_h=0
		#
		solturas=[8]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[750],[875]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[[0.9,0.9,0.3,0]]#taxa de qnts parasitoides saem do patch baixa e alta. 
		porcentagens=[0.15] #ocupacao inicial do host
		for solt_posi,solt_valor in enumerate(solturas):
			for pop_host_ini_porcentagem in porcentagens:
				porcentagem_popinihost=[pop_host_ini_porcentagem,0]
				for fracao_indv_migrante in indv_migrante:
					for borda_meio,hec10_BordaMeio in enumerate(hec10[solt_posi]): #pega borda e meio
						
						tempo=sorteio_pop_inicial_host(lin,col,porcentagem_popinihost,inicial_h,inicial_n,solt_valor,loop,lista_patches)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]
						
						tempo=alocar_parasitoide(borda_meio,solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,hec10_bordas,inicial_n,inicial_h,para_lista_inicial,loop)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]	
						patches_sorteados=tempo[2]				
						if borda_meio==0:
							bOUm="borda"
						else:
							bOUm="meio"
						wd_novo=f'{wd_original}\\{bOUm}'
						os.makedirs(wd_novo)
						os.chdir(wd_novo)
						zzzz=0
						replicas=0
						total_replicas=0
						while replicas<10 or total_replicas<50:
							modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,patches_sorteados,bOUm,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,borda_meio,loop,lib_extras,Pnumerico,Q1numerico,Q2numerico)

							ocupacao_file=open(f'{zzzz}_ocupacao_total.txt',"rb")
							zzzz+=1
							ocupacao_pck=ocupacao_file.read()
							ocupacao_file.close()
							ocupacao=pck.loads(ocupacao_pck)
							total_replicas+=1
							if total_replicas>49:
								break
							if ocupacao[0][30]>0 and ocupacao[1][30]>0:
								replicas+=1
	elif loop==26: #só Cotesia no dois pontos ao mesmo tempo com 6k cotesia total
		Pnumerico=52.5
		Q1numerico=30
		Q2numerico=0
		replicas=0 #replicas
		TS=[4,4,0] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,0]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=0
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=0
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.25,0,0.05,0]#[0.125,0.125,0.05,0]  #PQNH
		tempo_final=31
		inicial_n=41
		inicial_h=0
		#
		solturas=[8]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[375],[0]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[[0.9,0.9,0.3,0]]#taxa de qnts parasitoides saem do patch baixa e alta. 
		porcentagens=[0.15] #ocupacao inicial do host
		for solt_posi,solt_valor in enumerate(solturas):
			for pop_host_ini_porcentagem in porcentagens:
				porcentagem_popinihost=[pop_host_ini_porcentagem,0]
				for fracao_indv_migrante in indv_migrante:
					for borda_meio,hec10_BordaMeio in enumerate(hec10[solt_posi]): #pega borda e meio
						
						tempo=sorteio_pop_inicial_host(lin,col,porcentagem_popinihost,inicial_h,inicial_n,solt_valor,loop,lista_patches)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]
						
						tempo=alocar_parasitoide(borda_meio,solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,hec10_bordas,inicial_n,inicial_h,para_lista_inicial,loop)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]	
						patches_sorteados=tempo[2]				
						if borda_meio==0:
							bOUm="borda"
						else:
							bOUm="meio"
						wd_novo=f'{wd_original}\\{bOUm}'
						os.makedirs(wd_novo)
						os.chdir(wd_novo)
						zzzz=0
						total_replicas=0
						replicas=0
						while replicas<10 or total_replicas<50:
							modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,patches_sorteados,bOUm,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,borda_meio,loop,lib_extras,Pnumerico,Q1numerico,Q2numerico)

							ocupacao_file=open(f'{zzzz}_ocupacao_total.txt',"rb")
							zzzz+=1
							ocupacao_pck=ocupacao_file.read()
							ocupacao_file.close()
							ocupacao=pck.loads(ocupacao_pck)
							total_replicas+=1
							if total_replicas>49:
								break
							if ocupacao[0][30]>0:
								replicas+=1	
	elif loop==27: #só Cotesia no dois pontos ao mesmo tempo com 13k cotesia total
		Pnumerico=52.5
		Q1numerico=30
		Q2numerico=0
		replicas=0 #replicas
		TS=[4,4,0] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,0]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=0
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=0
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.25,0,0.05,0]#[0.125,0.125,0.05,0]  #PQNH
		tempo_final=31
		inicial_n=41
		inicial_h=0
		#
		solturas=[8]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[812.5],[0]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[[0.9,0.9,0.3,0]]#taxa de qnts parasitoides saem do patch baixa e alta. 
		porcentagens=[0.15] #ocupacao inicial do host
		for solt_posi,solt_valor in enumerate(solturas):
			for pop_host_ini_porcentagem in porcentagens:
				porcentagem_popinihost=[pop_host_ini_porcentagem,0]
				for fracao_indv_migrante in indv_migrante:
					for borda_meio,hec10_BordaMeio in enumerate(hec10[solt_posi]): #pega borda e meio
						
						tempo=sorteio_pop_inicial_host(lin,col,porcentagem_popinihost,inicial_h,inicial_n,solt_valor,loop,lista_patches)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]
						
						tempo=alocar_parasitoide(borda_meio,solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,hec10_bordas,inicial_n,inicial_h,para_lista_inicial,loop)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]	
						patches_sorteados=tempo[2]				
						if borda_meio==0:
							bOUm="borda"
						else:
							bOUm="meio"
						wd_novo=f'{wd_original}\\{bOUm}'
						os.makedirs(wd_novo)
						os.chdir(wd_novo)
						zzzz=0
						total_replicas=0
						replicas=0
						while replicas<10 or total_replicas<50:
							modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,patches_sorteados,bOUm,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,borda_meio,loop,lib_extras,Pnumerico,Q1numerico,Q2numerico)

							ocupacao_file=open(f'{zzzz}_ocupacao_total.txt',"rb")
							zzzz+=1
							ocupacao_pck=ocupacao_file.read()
							ocupacao_file.close()
							ocupacao=pck.loads(ocupacao_pck)
							total_replicas+=1
							if total_replicas>49:
								break
							if ocupacao[0][30]>0:
								replicas+=1	
	elif loop==28: #só Tetrastichus no dois pontos ao mesmo tempo com 13k tetrastichus total
		Pnumerico=52.5
		Q1numerico=30
		Q2numerico=0
		replicas=0 #replicas
		TS=[4,4,0] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,0]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=0
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=0
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0,0.25,0.05,0]#[0.125,0.125,0.05,0]  #PQNH
		tempo_final=31
		inicial_n=41
		inicial_h=0
		#
		solturas=[8]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[0],[812.5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[[0.9,0.9,0.3,0]]#taxa de qnts parasitoides saem do patch baixa e alta. 
		porcentagens=[0.15] #ocupacao inicial do host
		for solt_posi,solt_valor in enumerate(solturas):
			for pop_host_ini_porcentagem in porcentagens:
				porcentagem_popinihost=[pop_host_ini_porcentagem,0]
				for fracao_indv_migrante in indv_migrante:
					for borda_meio,hec10_BordaMeio in enumerate(hec10[solt_posi]): #pega borda e meio
						
						tempo=sorteio_pop_inicial_host(lin,col,porcentagem_popinihost,inicial_h,inicial_n,solt_valor,loop,lista_patches)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]
						
						tempo=alocar_parasitoide(borda_meio,solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,hec10_bordas,inicial_n,inicial_h,para_lista_inicial,loop)
						patch_iniciais=tempo[0]
						pop_iniciais=tempo[1]	
						patches_sorteados=tempo[2]				
						if borda_meio==0:
							bOUm="borda"
						else:
							bOUm="meio"
						wd_novo=f'{wd_original}\\{bOUm}'
						os.makedirs(wd_novo)
						os.chdir(wd_novo)
						zzzz=0
						total_replicas=0
						replicas=0
						while replicas<10 or total_replicas<50:
							modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,fracao_indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,patches_sorteados,bOUm,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,borda_meio,loop,lib_extras,Pnumerico,Q1numerico,Q2numerico)

							ocupacao_file=open(f'{zzzz}_ocupacao_total.txt',"rb")
							zzzz+=1
							ocupacao_pck=ocupacao_file.read()
							ocupacao_file.close()
							ocupacao=pck.loads(ocupacao_pck)
							total_replicas+=1
							if total_replicas>49:
								break
							if ocupacao[0][30]>0:
								replicas+=1	
	return
############################
############################
############################





lib_extras=False


loop=23
lin=50
col=50
distancia_maxima=3
lista_distancia=checando_existencia_lista_distancias(lin,col,distancia_maxima)
tempoo=checando_existencia_lista_viz_e_lista_de_patches(lin,col,distancia_maxima)
viz=tempoo[0]
lista_patches=tempoo[1]
lista_borda=tempoo[2]
if loop==23:
	wd_pasta_origem=f'{os.getcwd()}'
	for loop_numero in [loop,loop+1,loop+2,loop+3,loop+4,loop+5]:
		if loop_numero==23:
			pasta=f'{wd_pasta_origem}\\so cotesia'
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==24:
			#pasta=f'{wd_pasta_origem}\\ambos juntos'
			#os.makedirs(pasta)
			#os.chdir(pasta)
			continue
		if loop_numero==25:
			pasta=f'{wd_pasta_origem}\\ambos separados'
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==26:
			pasta=f'{wd_pasta_origem}\\cotesia nos dois6k'
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==27:
			pasta=f'{wd_pasta_origem}\\cotesia nos dois13k'
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==28:
			pasta=f'{wd_pasta_origem}\\tetrastichus nos dois13k'
			os.makedirs(pasta)
			os.chdir(pasta)
		loops_cenarios(loop_numero,viz,lista_patches,lista_borda)