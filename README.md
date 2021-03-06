#`Introduction`

The `AnimalGene2QTL` provides an interface to retrieve QTL data ,
gene data , SNP data.The package enables retrieval of large amounts
of data in a uniform way without the need to know the underlying
database schemas or write complex SQL queries.

#`Quick start`

## `Viewing AnimalGene2QTL dataset`
Every analysis with `AnimalGene2QTL` starts with Viewing a
`AnimalGene2QTL` dataset to use.
The function listQTL will display all available QTL dataset.

```{r eval=TRUE}
library(AnimalGene2QTL)
listQTL()
```

## `How to build a QTL query`
The `getAnimalQTL` function has six arguments that need to be
introduced:
*gene_filters*, *qtl_attributes*, *gene_values*, *data_set*
, *snp*, *snp_attributes*.

* **gene_filters**: define a restriction on the query.
For example you want to restrict the output to all QTL located
on the human X chromosome
then the `gene_filters`:`chromosome_name` can be used with gene
values `X`.
The `listGeneAF` function displays all available gene filters
in the selected dataset.
* **qtl_attributes**: define the values we are interested in to retrieve.
For example we want to retrieve the QTL ID or QTL name.
The `listQTLAF` function displays all available QTL attributes
in the selected dataset.
* **gene_values**: The input gene data.
* **data_set**: There are five numbers:1,2,3,4,5.`listQTL`function display
which to input.
* **snp**: default FALSE,when it is TRUE,that mean retrieve SNP.
* **snp_attributes**: input the attribute of SNP what to retrieve,`listSNPattributes`function display it.
```{r eval=TRUE}
attributes <- listQTLAF()
head(attributes)
filters <- listGeneAF(1)
head(filters)
```

### `getAnimalQTL() function`
The `getAnimalQTL` function is the main query function in
`AnimalGene2QTL`.
It has four main arguments.

* **qtl_attributes**: is a vector of attributes that one
wants to retrieve (= the output of the query).
* **gene_filters**: is a vector of filter that one wil use
as input to the query.
* **gene_values**: a vector of values for the *gene_filters*.
In case multple filter are in use, the *gene_values* argument
requires a list of values where each position
in the list corresponds to the position of the filters in the
*gene_filters* argument (see examples below).
* **data_set**: choose which animal QTL you want to retrieve.

Now that we selected a `AnimalGene2QTL` dataset, and know
about `qtl_attributes`,
`gene_filters`, and the `gene_values` for `gene_filters`;
we can build a `AnimalGene2QTL`
query. Let us make an easy query for the following problem:
We have a list of gene identifiers from the Ensembl and we
want to retrieve the QTL identifiers.
Let us now run the query:

```{r eval=TRUE}
geneid <- c("ENSBTAG00000009851", "ENSBTAG00000005101", "ENSBTAG00000036262")
qtl <- getAnimalQTL(qtl_attributes=c('QTL_ID'),gene_filters='ensembl_gene_id'
,gene_values=geneid,data_set=1);
head(qtl);
```

## `Examples of AnimalGene2QTL queries`
In the sections below a variety of example
queries are described.
Every example is written
as a task, and we have to come up with a `AnimalGene2QTL`
solution to the problem.

### `Task 1: Retrieve all identifiers of genes by QTL`
`identifiers:"64577","2199","2354"`

```{r eval=TRUE}
qtlid <- c("64577","2199","2354");
gene <- getAnimalGene(gene_attributes="ensembl_gene_id",qtl_filters=c("QTL_ID")
,qtl_values=qtlid, data_set=2);
head(gene)
```

### `Task 2: Retrieve all identifiers of SNP by QTL`
`identifiers:"4097"`

```{r eval=TRUE}
snp <- getSNPbyQTL('refsnp_id','QTL_ID','4097',2);
head(snp)
```