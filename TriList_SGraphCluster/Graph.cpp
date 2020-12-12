#include "Utility.h"
#include "Graph.h"

Graph::Graph(const char *_dir) {
	dir = string(_dir);

	n = m = 0;

	pstart = NULL;
	edges = NULL;
}

Graph::~Graph() {
	if(pstart != NULL) {
		delete[] pstart;
		pstart = NULL;
	}
	if(edges != NULL) {
		delete[] edges;
		edges = NULL;
	}
}

void Graph::read_graph() {
	FILE *f = open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != (int)sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	// printf("\tn = %u; m = %u\n", n, m/2);

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, f);

#ifdef _DEBUG_
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	if(sum != m) printf("WA input graph\n");
#endif

	fclose(f);

	f = open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart == NULL) pstart = new ui[n+1];
	if(edges == NULL) edges = new int[m];

	ui *buf = new ui[n];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		//printf("%d %d\n", i, degree[i]);
		if(degree[i] > 0) fread(buf, sizeof(int), degree[i], f);

		for(ui j = 0;j < degree[i];j ++) edges[pstart[i] + j] = buf[j];

		pstart[i+1] = pstart[i] + degree[i];
	}

	delete[] buf;
	buf = NULL;

	fclose(f);

#ifdef _DEBUG_
	for(ui i = 0;i < n;i ++) {
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(edges[j] == i) {
				printf("Self loop\n");
				//exit(1);
			}

			if(j > pstart[i]&&edges[j] <= edges[j-1]) {
				printf("Edges not sorted in increasing id order!\n");
				//exit(1);
			}
		}
	}
#endif

	delete[] degree;
}

ui Graph::binary_search(const int *array, ui b, ui e, int val) {
#ifdef _DEBUG_
	if(e <= b) printf("??? WA1 in binary_search\n");
#endif
	-- e;
	if(array[e] < val) return e+1;
	while(b < e) {
		ui mid = b + (e-b)/2;
		if(array[mid] >= val) e = mid;
		else b = mid+1;
	}
#ifdef _DEBUG_
	if(array[e] < val) printf("??? WA2 in binary_search\n");
#endif
	return e;
}

void Graph::scan(const char *eps_s, const char *mu_s, int out) {
#ifdef _LINUX_
	struct timeval start;
	gettimeofday(&start, NULL);
#else
	int start = clock();
#endif

	int eps_a2 = 0, eps_b2 = 0;
	get_eps(eps_s, eps_a2, eps_b2);
	eps_a2 *= eps_a2, eps_b2 *= eps_b2;

	int mu = atoi(mu_s);

	ui *reverse = new ui[m];
	cross_link(reverse);

#ifdef _LINUX_
	struct timeval end1;
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#else
	int end1 = clock();
#endif

	ui *deg = new ui[n];
	for(ui i = 0;i < n;i ++) deg[i] = pstart[i+1]-pstart[i];

	ui *adj = new ui[n];
	memset(adj, 0, sizeof(ui)*n);

	ui *similar = new ui[m];
	memset(similar, 0, sizeof(ui)*m);

	ui *vv = new ui[n];

	for(ui u = 0;u < n;u ++) {
		ui vv_n = 0;
		for(ui j = pstart[u];j < pstart[u+1];j ++) {
			ui v = edges[j];
			if(deg[v] < deg[u]||(deg[v] == deg[u]&&v < u)) {
				adj[v] = j+1;
				vv[vv_n ++] = j;
			}
		}

		for(ui j = 0;j < vv_n;j ++) {
			ui v = edges[vv[j]];

			for(ui k = pstart[v];k < pstart[v+1]&&edges[k] < v;k ++) if(adj[edges[k]]) {
				ui pos = vv[j];
				++ similar[pos]; ++ similar[reverse[pos]];
				++ similar[k]; ++ similar[reverse[k]];
				pos = adj[edges[k]] - 1;
				++ similar[pos]; ++ similar[reverse[pos]];
			}
		}

		for(ui j = 0;j < vv_n;j ++) adj[edges[vv[j]]] = 0;
	}

	for(ui u = 0;u < n;u ++) {
		for(ui j = pstart[u];j < pstart[u+1];j ++) {
			int v = edges[j];
			if(v < u) continue;

			similar[j] += 2;

			if(((long long)similar[j])*((long long)similar[j])*eps_b2 >= ((long long)(deg[u]+1))*((long long)(deg[v]+1))*eps_a2) similar[j] = 1;
			else similar[j] = 0;

			similar[reverse[j]] = similar[j];
		}
	}

	delete[] vv; vv = NULL;
	delete[] deg; deg = NULL;
	delete[] adj; adj = NULL;

	int *similar_degree = new int[n];
	memset(similar_degree, 0, sizeof(int)*n);

	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		if(similar[j] == 1) ++ similar_degree[i];
	}

	int *pa = new int[n];
	int *rank = new int[n];
	for(ui i = 0;i < n;i ++) {
		pa[i] = i;
		rank[i] = 0;
	}

	for(ui i = 0;i < n;i ++) {
		if(similar_degree[i] < mu) continue;
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(similar_degree[edges[j]] < mu) continue;

			if(similar[j] == 1) my_union(pa, rank, i, edges[j]);
		}
	}

	int *id = new int[n];
	for(ui i = 0;i < n;i ++) id[i] = n;

	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		int x = find_root(pa, i);
		if(i < id[x]) id[x] = i;
	}

	vector<pair<int,int> > noncore_cluster;
	noncore_cluster.reserve(n);

	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
		for(ui j = pstart[i];j < pstart[i+1];j ++) {
			if(similar_degree[edges[j]] >= mu) continue;

			if(similar[j] == 1) noncore_cluster.pb(mp(id[pa[i]], edges[j]));
		}
	}

#ifdef _LINUX_
	struct timeval end;
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - end1.tv_sec;
	useconds = end.tv_usec - end1.tv_usec;
	mtime = seconds*1000000 + useconds;

	//printf("Preprocess time: %lld\nCluster time: %lld\n", mtime1, mtime);
#else
	int end = clock();

	printf("Preprocess time: %d\nCluster time: %d\n", end1-start,end-end1);
#endif

	if(out) output(similar_degree, id, pa, noncore_cluster, eps_s, mu_s);

	delete[] id;
	delete[] pa;
	delete[] rank;
	delete[] similar;
	delete[] similar_degree;
	delete[] reverse;
}

void Graph::output(int *similar_degree, int *id, int *pa, vector<pair<int,int> > &noncore_cluster, const char *eps_s, const char *mu) {
	printf("\tStart write result into disk\n");

	string out_name = dir + "/result-" + string(eps_s) + "-" + string(mu) + ".txt";
	FILE *fout = open_file(out_name.c_str(), "w");

	int miu = atoi(mu);
	for(ui i = 0;i < n;i ++) if(similar_degree[i] >= miu) {
		fprintf(fout, "c %d %d\n", i, id[find_root(pa, i)]);
	}

	sort(noncore_cluster.begin(), noncore_cluster.end());
	noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
	for(ui i = 0;i < noncore_cluster.size();i ++) {
		fprintf(fout, "n %d %d\n", noncore_cluster[i].second, noncore_cluster[i].first);
	}

	fclose(fout);
}

void Graph::my_union(int *pa, int *rank, int u, int v) {
	int ru = find_root(pa, u);
	int rv = find_root(pa, v);

	if(ru == rv) return ;

	if(rank[ru] < rank[rv]) pa[ru] = rv;
	else if(rank[ru] > rank[rv]) pa[rv] = ru;
	else {
		pa[rv] = ru;
		++ rank[ru];
	}
}

int Graph::find_root(int *pa, int u) {
	int x = u;
	while(pa[x] != x) x = pa[x];

	while(pa[u] != x) {
		int tmp = pa[u];
		pa[u] = x;
		u = tmp;
	}

	return x;
}

void Graph::cross_link(ui *reverse) {
	for(ui i = 0;i < n;i ++) for(ui j = pstart[i];j < pstart[i+1];j ++) {
		int v = edges[j];
		int d1 = pstart[i+1] - pstart[i], d2 = pstart[v+1] - pstart[v];

		if(d1 < d2||(d1 == d2&&i < v)) continue;

		ui r_id = binary_search(edges, pstart[v], pstart[v+1], i);
#ifdef _DEBUG_
		if(r_id >= pstart[v+1]||i != edges[r_id]) printf("??? WA in cross_link\n");
#endif

		reverse[j] = r_id;
		reverse[r_id] = j;
	}
}

void Graph::get_eps(const char *eps_s, int &eps_a, int &eps_b) {
	int i = 0;
	eps_a = 0; eps_b = 1;
	while(eps_s[i] != '\0'&&eps_s[i] != '.') {
		eps_a = eps_a*10 + (eps_s[i]-'0');
		++ i;
	}

	if(eps_s[i] == '.') {
		++ i;
		while(eps_s[i] != '\0') {
			eps_a = eps_a*10 + (eps_s[i]-'0');
			eps_b *= 10;
			++ i;
		}
	}

	if(eps_a > eps_b||eps_b > 100||eps_a <= 0) {
		printf("??? Wrong eps format: %d/%d, %s\n", eps_a, eps_b, eps_s);
		exit(1);
	}
}
