import multiprocessing
import os
import matplotlib
import umap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shapely.geometry as geom

from statsmodels.nonparametric.smoothers_lowess import lowess
from pandas.api.types import is_string_dtype
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.manifold import LocallyLinearEmbedding, TSNE, SpectralEmbedding

def project_point_to_curve_distance(XP,p):
    curve = geom.LineString(XP)
    point = geom.Point(p)
    #distance from point to curve
    dist_p_to_c = point.distance(curve)
    return dist_p_to_c

def select_variable_genes(
        adata,
        loess_frac=0.01,
        percentile=95,
        n_genes = None,
        n_jobs = multiprocessing.cpu_count(),
        save_fig=False,
        fig_name='std_vs_means.pdf',
        fig_path=None,
        fig_size=(4,4),
        pad=1.08,
        w_pad=None,
        h_pad=None
    ):

    """Select the most variable genes.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    loess_frac: `float`, optional (default: 0.1)
        Between 0 and 1. The fraction of the data used when estimating each y-value in LOWESS function.
    percentile: `int`, optional (default: 95)
        Between 0 and 100. Specify the percentile to select genes.Genes are ordered based on its distance from the fitted curve.
    n_genes: `int`, optional (default: None)
        Specify the number of selected genes. Genes are ordered based on its distance from the fitted curve.
    n_jobs: `int`, optional (default: all available cpus)
        The number of parallel jobs to run when calculating the distance from each gene to the fitted curve
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_size: `tuple`, optional (default: (4,4))
        figure size.
    fig_path: `str`, optional (default: '')
        if empty, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'std_vs_means.pdf')
        if save_fig is True, specify figure name.
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.

    Returns
    -------
    updates `adata` with the following fields.
    var_genes: `numpy.ndarray` (`adata.obsm['var_genes']`)
        Store #observations × #var_genes data matrix used for subsequent dimension reduction.
    var_genes: `pandas.core.indexes.base.Index` (`adata.uns['var_genes']`)
        The selected variable gene names.
    """

    # if(not issparse(adata.X)):
    #     adata.X = csr_matrix(adata.X) 

    if(fig_path is None):
        fig_path = adata.uns['workdir']  
    fig_size = matplotlib.rcParams['figure.figsize'] if fig_size is None else fig_size
    mean_genes = np.mean(adata.X,axis=0)
    std_genes = np.std(adata.X,ddof=1,axis=0)
    loess_fitted = lowess(std_genes,mean_genes,return_sorted=False,frac=loess_frac)
    residuals = std_genes - loess_fitted
    XP = np.column_stack((np.sort(mean_genes),loess_fitted[np.argsort(mean_genes)]))
    mat_p = np.column_stack((mean_genes,std_genes))
    with multiprocessing.Pool(processes=n_jobs) as pool:
        dist_point_to_curve = pool.starmap(project_point_to_curve_distance,[(XP,mat_p[i,]) for i in range(XP.shape[0])])
    mat_sign = np.ones(XP.shape[0])
    mat_sign[np.where(residuals<0)[0]] = -1
    dist_point_to_curve = np.array(dist_point_to_curve)*mat_sign
    if(n_genes is None):
        cutoff = np.percentile(dist_point_to_curve,percentile)
        id_var_genes = np.where(dist_point_to_curve>cutoff)[0]
        id_non_var_genes = np.where(residuals<=cutoff)[0]
    else:
        id_var_genes = np.argsort(dist_point_to_curve)[::-1][:n_genes]
        id_non_var_genes = np.argsort(dist_point_to_curve)[::-1][n_genes:]

    adata.obsm['var_genes'] = adata.X[:,id_var_genes].copy()
    adata.uns['var_genes'] = adata.var_names[id_var_genes]
    print(str(len(id_var_genes))+' variable genes are selected')
    ###plotting
    fig = plt.figure(figsize=fig_size)      
    plt.scatter(mean_genes[id_non_var_genes], std_genes[id_non_var_genes],s=5,alpha=0.2,zorder=1,c='#6baed6')
    plt.scatter(mean_genes[id_var_genes], std_genes[id_var_genes],s=5,alpha=0.9,zorder=2,c='#EC4E4E')
    plt.plot(np.sort(mean_genes), loess_fitted[np.argsort(mean_genes)],linewidth=3,zorder=3,c='#3182bd')
    plt.xlabel('mean value')
    plt.ylabel('standard deviation')
    plt.locator_params(axis='x',nbins=5)
    plt.locator_params(axis='y',nbins=5)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
    if(save_fig):
        plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
        plt.close(fig)
    return None


def plot_visualization_2D(adata,method='umap',n_neighbors=50, nb_pct=None,perplexity=30.0,color=None,use_precomputed=True,
                          fig_size=None,fig_ncol=3,fig_legend_ncol=1,fig_legend_order = None,
                          vmin=None,vmax=None,alpha=0.8,
                          pad=1.08,w_pad=None,h_pad=None,
                          save_fig=False,fig_path=None,fig_name='visualization_2D.pdf',
                          plotly=False):

    """ Visualize the results in 2D plane
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix. 
    method: `str`, optional (default: 'umap')
        Choose from {{'umap','tsne'}}
        Method used for visualization.
        'umap': Uniform Manifold Approximation and Projection      
        'tsne': t-Distributed Stochastic Neighbor Embedding
    n_neighbors: `int`, optional (default: 50)
        The number of neighbor cells (only valid when 'umap' is specified).
    nb_pct: `float`, optional (default: None)
        The percentage of neighbor cells (when sepcified, it will overwrite n_neighbors).
    perplexity: `float`, optional (default: 30.0)
        The perplexity used for tSNE. 
    color: `list` optional (default: None)
        Column names of observations (adata.obs.columns) or variable names(adata.var_names). A list of names to be plotted.    
    use_precomputed: `bool`, optional (default: True)
        If True, the visualization coordinates from previous computing will be used
    fig_size: `tuple`, optional (default: None)
        figure size.
    fig_legend_order: `dict`,optional (default: None)
        Specified order for the appearance of the annotation keys.Only valid for ategorical variable  
        e.g. fig_legend_order = {'ann1':['a','b','c'],'ann2':['aa','bb','cc']}
    fig_legend_ncol: `int`, optional (default: 1)
        The number of columns that the legend has.
    vmin,vmax: `float`, optional (default: None)
        The min and max values are used to normalize continuous values. If None, the respective min and max of continuous values is used.
    alpha: `float`, optional (default: 0.8)
        0.0 transparent through 1.0 opaque
    pad: `float`, optional (default: 1.08)
        Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    h_pad, w_pad: `float`, optional (default: None)
        Padding (height/width) between edges of adjacent subplots, as a fraction of the font size. Defaults to pad.
    save_fig: `bool`, optional (default: False)
        if True,save the figure.
    fig_path: `str`, optional (default: None)
        if save_fig is True, specify figure path. if None, adata.uns['workdir'] will be used.
    fig_name: `str`, optional (default: 'visualization_2D.pdf')
        if save_fig is True, specify figure name.
    plotly: `bool`, optional (default: False)
        if True, plotly will be used to make interactive plots

    Returns
    -------
    updates `adata` with the following fields. (Depending on `method`)
    vis_trans_umap : `umap.UMAP` (`adata.uns['vis_trans_umap']`)
        Store umap object
    vis_trans_tsne : `sklearn.manifold._t_sne.TSNE` (`adata.uns['vis_trans_tsne']`)
        Store tsne object
    X_vis_umap: `numpy.ndarray` (`adata.obsm['X_vis_umap']`)
        Store #observations × 2 umap data matrix. 
    X_vis_tsne: `numpy.ndarray` (`adata.obsm['X_vis_tsne']`)
        Store #observations × 2 tsne data matrix. 
    X_vis : `numpy.ndarray` (`adata.obsm['X_vis']`)
        A #observations × 2 data matrix after visualization.   
    """    
    if(fig_path is None):
        fig_path = adata.uns['workdir']
    fig_size = matplotlib.rcParams['figure.figsize'] if fig_size is None else fig_size

    if(method not in ['umap','tsne']):
        raise ValueError("unrecognized method '%s'" % method)
    if(color is None):
        color = ['label']
    ###remove duplicate keys
    color = list(dict.fromkeys(color))

    dict_ann = dict()
    for ann in color:
        if(ann in adata.obs.columns):
            dict_ann[ann] = adata.obs[ann]
        elif(ann in adata.var_names):
            dict_ann[ann] = adata.obs_vector(ann)
        else:
            raise ValueError("could not find '%s' in `adata.obs.columns` and `adata.var_names`"  % (ann))
    input_data = adata.obsm['X_dr']
    if(nb_pct!=None):
        n_neighbors = int(np.around(input_data.shape[0]*nb_pct)) 
    if(method == 'umap'):       
        if(use_precomputed and ('X_vis_umap' in adata.obsm_keys())):
            print('Importing precomputed umap visualization ...')
            embedding = adata.obsm['X_vis_umap']
        else:
            reducer = umap.UMAP(n_neighbors=n_neighbors,n_components=2,random_state=42)
            trans = reducer.fit(input_data)
            embedding = trans.embedding_
            adata.uns['vis_trans_umap'] = trans
            adata.obsm['X_vis_umap'] = embedding
    if(method == 'tsne'):
        if(use_precomputed and ('X_vis_tsne' in adata.obsm_keys())):
            print('Importing precomputed tsne visualization ...')
            embedding = adata.obsm['X_vis_tsne']
        else:
            reducer = TSNE(n_components=2, init='pca',perplexity=perplexity, random_state=42)
            trans = reducer.fit(input_data)
            embedding = trans.embedding_
            adata.uns['vis_trans_tsne'] = trans
            adata.obsm['X_vis_tsne'] = embedding
    adata.obsm['X_vis'] = embedding
    if('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['plot_visualization_2D'] = {'method':method}
    
    df_plot = pd.DataFrame(index=adata.obs.index,data = embedding,columns=[method.upper()+str(x) for x in [1,2]])
    for ann in color:
        df_plot[ann] = dict_ann[ann]
    df_plot_shuf = df_plot.sample(frac=1,random_state=100)
    
    legend_order = {ann:np.unique(df_plot_shuf[ann]) for ann in color if is_string_dtype(df_plot_shuf[ann])}
    if(fig_legend_order is not None):
        if(not isinstance(fig_legend_order, dict)):
            raise TypeError("`fig_legend_order` must be a dictionary")
        for ann in fig_legend_order.keys():
            if(ann in legend_order.keys()):
                legend_order[ann] = fig_legend_order[ann]
            else:
                print("'%s' is ignored for ordering legend labels due to incorrect name or data type" % ann)

    # if(plotly):
    #     for ann in color:
    #         fig = px.scatter(df_plot_shuf, x='Dim'+str(comp1+1), y='Dim'+str(comp2+1),color=ann,
    #                             opacity=alpha,width=500,height=500,
    #                             color_continuous_scale=px.colors.sequential.Viridis,
    #                             color_discrete_map=adata.uns[ann+'_color'] if ann+'_color' in adata.uns_keys() else {})
    #         fig.update_layout(legend= {'itemsizing': 'constant'}) 
    #         fig.show(renderer="notebook")
    else:
        if(len(color)<fig_ncol):
            fig_ncol=len(color)
        fig_nrow = int(np.ceil(len(color)/fig_ncol))
        fig = plt.figure(figsize=(fig_size[0]*fig_ncol*1.05,fig_size[1]*fig_nrow))
        for i,ann in enumerate(color):
            ax_i = fig.add_subplot(fig_nrow,fig_ncol,i+1)
            if(is_string_dtype(df_plot[ann])):
                sc_i=sns.scatterplot(ax=ax_i,
                                    x=method.upper()+'1', y=method.upper()+'2', 
                                    hue=ann,hue_order = legend_order[ann],
                                    data=df_plot_shuf,
                                    alpha=alpha,linewidth=0,
                                    palette= adata.uns[ann+'_color'] \
                                            if (ann+'_color' in adata.uns_keys()) and (set(adata.uns[ann+'_color'].keys()) >= set(np.unique(df_plot_shuf[ann]))) \
                                            else None
                                    )
                legend_handles, legend_labels = ax_i.get_legend_handles_labels()
                ax_i.legend(handles=legend_handles, labels=legend_labels,
                            bbox_to_anchor=(1, 0.5), loc='center left', ncol=fig_legend_ncol,
                            frameon=False,
                            )
                if(ann+'_color' not in adata.uns_keys()):
                    colors_sns = sc_i.get_children()[0].get_facecolors()
                    colors_sns_scaled = (255*colors_sns).astype(int)
                    adata.uns[ann+'_color'] = {df_plot_shuf[ann][i]:'#%02x%02x%02x' % (colors_sns_scaled[i][0], colors_sns_scaled[i][1], colors_sns_scaled[i][2])
                                                for i in np.unique(df_plot_shuf[ann],return_index=True)[1]}
                ### remove legend title
                # ax_i.get_legend().texts[0].set_text("")
            else:
                vmin_i = df_plot[ann].min() if vmin is None else vmin
                vmax_i = df_plot[ann].max() if vmax is None else vmax
                sc_i = ax_i.scatter(df_plot_shuf[method.upper()+'1'], df_plot_shuf[method.upper()+'2'],
                                    c=df_plot_shuf[ann],vmin=vmin_i,vmax=vmax_i,alpha=alpha)
                cbar = plt.colorbar(sc_i,ax=ax_i, pad=0.01, fraction=0.05, aspect=40)
                cbar.solids.set_edgecolor("face")
                cbar.ax.locator_params(nbins=5)                    
            ax_i.set_xlabel(method.upper()+'1')
            ax_i.set_ylabel(method.upper()+'2',labelpad=2)
            ax_i.get_xaxis().set_ticks([])
            ax_i.get_yaxis().set_ticks([])
            ax_i.set_title(ann)
#             plt.subplots_adjust(hspace=hspace,wspace=wspace)
        plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)
        if(save_fig):
            plt.savefig(os.path.join(fig_path,fig_name),pad_inches=1,bbox_inches='tight')
            plt.close(fig)

def dimension_reduction(adata,n_neighbors=50, nb_pct = None,n_components = 3,n_jobs = 1,
                        feature='var_genes',method = 'se',eigen_solver=None):

    """Perform dimension reduction.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    n_neighbors: `int`, optional (default: 50)
        The number of neighbor cells used for manifold learning (only valid when 'mlle','se', or 'umap' is specified).
    nb_pct: `float`, optional (default: None)
        The percentage of neighbor cells (when sepcified, it will overwrite n_neighbors).
    n_components: `int`, optional (default: 3)
        Number of components to keep.
    n_jobs: `int`, optional (default: 1)
        The number of parallel jobs to run.
    feature: `str`, optional (default: 'var_genes')
        Choose from {{'var_genes','top_pcs','all'}}
        Feature used for dimension reduction.
        'var_genes': most variable genes
        'top_pcs': top principal components
        'all': all available features (genes)
    method: `str`, optional (default: 'se')
        Choose from {{'se','mlle','umap','pca'}}
        Method used for dimension reduction.
        'se': Spectral embedding algorithm
        'mlle': Modified locally linear embedding algorithm
        'umap': Uniform Manifold Approximation and Projection
        'pca': Principal component analysis
    eigen_solver: `str`, optional (default: None)
        For 'mlle', choose from {{'arpack', 'dense'}}
        For 'se', choose from {{'arpack', 'lobpcg', or 'amg'}}
        The eigenvalue decomposition strategy to use
   
    Returns
    -------
    updates `adata` with the following fields.
    
    X_dr : `numpy.ndarray` (`adata.obsm['X_dr']`)
        A #observations × n_components data matrix after dimension reduction.
    X_mlle : `numpy.ndarray` (`adata.obsm['X_mlle']`)
        Store #observations × n_components data matrix after mlle.
    X_se : `numpy.ndarray` (`adata.obsm['X_se']`)
        Store #observations × n_components data matrix after spectral embedding.    
    X_umap : `numpy.ndarray` (`adata.obsm['X_umap']`)
        Store #observations × n_components data matrix after umap.
    X_pca : `numpy.ndarray` (`adata.obsm['X_pca']`)
        Store #observations × n_components data matrix after pca.
    trans_mlle : `sklearn.manifold.locally_linear.LocallyLinearEmbedding` (`adata.uns['trans_mlle']`)
        Store mlle object
    trans_se : `sklearn.manifold.spectral_embedding_.SpectralEmbedding` (`adata.uns['trans_se']`)
        Store se object
    trans_umap : `umap.UMAP` (`adata.uns['trans_umap']`)
        Store umap object
    trans_pca : `sklearn.decomposition.PCA` (`adata.uns['trans_pca']`)
        Store pca object 
    """

    if(feature not in ['var_genes','top_pcs','all']):
        raise ValueError("unrecognized feature '%s'" % feature)
    if(method not in ['mlle','se','umap','pca']):
        raise ValueError("unrecognized method '%s'" % method)
    if(feature == 'var_genes'):
        input_data = adata.obsm['var_genes']
    if(feature == 'top_pcs'):
        input_data = adata.obsm['top_pcs']
    if(feature == 'all'):
        input_data = adata.X
    print('feature ' + feature + ' is being used ...')
    print(str(n_jobs)+' cpus are being used ...')
    if(nb_pct!=None):
        n_neighbors = int(np.around(input_data.shape[0]*nb_pct))

    if(method == 'mlle'):
        np.random.seed(2)
        if(eigen_solver==None):
            if(input_data.shape[0]<=2000):
                reducer = LocallyLinearEmbedding(n_neighbors=n_neighbors, 
                                                     n_components=n_components,
                                                     n_jobs = n_jobs,
                                                     method = 'modified',eigen_solver = 'dense',random_state=10,
                                                     neighbors_algorithm = 'kd_tree')
            else:
                reducer = LocallyLinearEmbedding(n_neighbors=n_neighbors, 
                                                     n_components=n_components,
                                                     n_jobs = n_jobs,
                                                     method = 'modified',eigen_solver = 'arpack',random_state=10,
                                                     neighbors_algorithm = 'kd_tree')
                
        else:
            reducer = LocallyLinearEmbedding(n_neighbors=n_neighbors, 
                                                 n_components=n_components,
                                                 n_jobs = n_jobs,
                                                 method = 'modified',eigen_solver = eigen_solver,random_state=10,
                                                 neighbors_algorithm = 'kd_tree')        
        trans = reducer.fit(input_data)
        adata.uns['trans_mlle'] = trans
        adata.obsm['X_mlle'] = trans.embedding_
        adata.obsm['X_dr'] = trans.embedding_
    if(method == 'se'):
        np.random.seed(2)
        reducer = SpectralEmbedding(n_neighbors=n_neighbors, 
                                         n_components=n_components,
                                         n_jobs = n_jobs,
                                         eigen_solver = eigen_solver,random_state=10)
        trans = reducer.fit(input_data)
        adata.uns['trans_se'] = trans
        adata.obsm['X_se'] = trans.embedding_
        adata.obsm['X_dr'] = trans.embedding_
    if(method == 'umap'):
        reducer = umap.UMAP(n_neighbors=n_neighbors,n_components=n_components,random_state=42)
        trans = reducer.fit(input_data)
        adata.uns['trans_umap'] = trans
        adata.obsm['X_umap'] = trans.embedding_
        adata.obsm['X_dr'] = trans.embedding_
    if(method == 'pca'):
        reducer = sklearnPCA(n_components=n_components,svd_solver='arpack',random_state=42)
        trans = reducer.fit(input_data)
        adata.uns['trans_pca'] = trans
        adata.obsm['X_pca'] = trans.transform(input_data) 
        adata.obsm['X_dr'] = adata.obsm['X_pca']
    if('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['dimension_reduction'] = {'feature':feature,'method':method,'n_components':n_components}
    return None
