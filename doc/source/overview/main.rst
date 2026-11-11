Main configuration file
==============================
The main configuration file is used to define the global settings, notably the paths to the different sources of data. The parameters of the main configuration file are described in Table 2.

Table 2: Configuration parameters of the main configuration file

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Name
     - Type
     - Description
   * - name
     - STRING
     - Name of the data processing. This name is used in the path of the saved products in the following manner:

       ``<saveimgbasepath>/<name>/<YYYY-MM-DD>/<datasetname>/<prodname>/<outputname>``
   * - datapath
     - STRING
     - Base directory of the rainbow raw data. This field must have a trailing '/'. The raw data files of a scan can be found using the following file path:

       ``<datapath>/<scanname>/<YYYY-MM-DD>/<YYYYMMDDHHMMSS00datatype>.<ext>``
   * - configpath
     - STRING
     - Absolute base directory of the configuration files. This directory contains clutter maps, filter coefficients, antenna pattern, and the data processing configuration files. If any path within a pyrad config file is provided as a relative path, it will be converted to an absolute path, by prepending it with the path specified with "configpath". If not provided, the directory of the main config file will be used.
   * - locationConfigFile
     - STRING
     - File nameof the location configuration file. If only the filename (and not directory) is provided it will be prepended with configpath. Described in Section 3.2.
   * - productConfigFile
     - STRING
     - File name of the product configuration file.  If only the filename (and not directory) is provided it will be prepended with configpath. Described in Section 4.
   * - lastStateFile
     - STRING
     - File name (with full path) of the file containing the time of the last processed scan. Used in particular for real-time processing.
   * - imgformat
     - STRING
       STRARR
     - File format(s) of the images. The following formats are supported: eps, png, and jpg. If ``saveimg`` is set to 0, this field is not used.
   * - saveimgbasepath
     - STRING
     - Base directory for the images to save. The directory structure looks as follows:

       ``<saveimgbasepath>/<name>/<YYYY-MM-DD>/<datasetname>/<prodname>/<outputname>``

       If ``saveimg`` is set to 0, this field is not used.
   * - s3EndpointRead
     - STRING
     - OPTIONAL. Url to an S3 endpoint containing input radar data. The format must be https://endpoint.domain (e.g. https://eu-central-1.linodeobjects.com/), the https:// is not mandatory.
   * - s3BucketRead
     - STRING
     - OPTIONAL. Name of an S3 bucket containing input radar data. It has to be used together with ``s3EndpointRead`` and ``s3PathRead`` to be able to read input radar data from a bucket. The procedure will only work by setting the environment variables S3_KEY_READ and S3_SECRET_READ.
   * - s3PathRead
     - STRING
     - OPTIONAL. Path where to find input radar data in an S3 bucket. The data will be retrieved at url https://s3BucketRead.s3EndpointRead/s3PathRead/*
   * - rm_s3_file
     - INT
     - OPTIONAL. If set, input radar data files downloaded from an S3 bucket will be removed after reading.
   * - s3EndpointWrite
     - STRING
     - OPTIONAL. Url to an S3 endpoint where to store output data. The format must be https://endpoint.domain (e.g. https://eu-central-1.linodeobjects.com/), the https:// is not mandatory.
   * - s3BucketWrite
     - STRING
     - OPTIONAL. Name of an S3 bucket for output radar data. It has to be used together with ``s3EndpointWrite`` and ``s3PathWrite`` to save output data to a bucket. Requires S3_KEY_WRITE and S3_SECRET_WRITE environment variables.
   * - s3PathWrite
     - STRING
     - OPTIONAL. Path where to save output radar data in an S3 bucket. The data will be saved at url https://s3BucketRead.s3EndpointRead/s3PathRead/filename
   * - s3SplitExtensionWrite
     - INT
     - OPTIONAL. If set to 1, files with multiple extensions are stored separately in the S3 bucket. For example, myfolder/out.csv becomes myfolder/csv/out.csv. Default is 0 (False).
   * - s3Verify
     - INT
     - OPTIONAL. If set to 0, disables SSL verification for S3 operations. Default is 1.
   * - s3Certificates
     - INT
     - OPTIONAL. Path to SSL certificates, only needed if custom certificates are used.
   * - MySQLhost
     - STRING
     - REQUIRED. Hostname or IP of the MySQL server. Use ``localhost`` when the database runs on the same machine.
   * - MySQLdatabase
     - STRING
     - OPTIONAL. Name of the MySQL database to connect to. The database must already exist.
   * - MySQLuser
     - STRING
     - OPTIONAL. Username used to authenticate with the MySQL server. The user must have appropriate privileges on the selected database.
   * - MySQLport
     - INTEGER
     - OPTIONAL. Port number of the MySQL server. Default is ``3306``. Use a different value (e.g. ``3307``) if the server is exposed on a custom port.
   * - loadbasepath
     - STRING
     - OPTIONAL. Base path of saved data. Defaults to ``saveimgbasepath``.
   * - loadname
     - STRING
     - OPTIONAL. Name of the saved data processing. Used for saved volume loading. Defaults to ``name``.
   * - gecsxbasepath
     - STRING
     - OPTIONAL. Base path of saved visibility fields generated by the GECSX tool.
   * - gecsxname
     - STRING
     - OPTIONAL. Name of the GECSX dataset containing the visibility fields. Used for GECSX volume loading.
   * - dempath
     - STRING
     - OPTIONAL. Base directory of the Digital Elevation Model (DEM) files. Used to load radar visibility.
   * - smnpath
     - STRING
     - OPTIONAL. Base directory of SwissMetNet stations data. Used for comparison with rain gauges.
   * - disdropath
     - STRING
     - OPTIONAL. Base directory of disdrometer data. Used for comparison with disdrometers.
   * - solarfluxpath
     - STRING
     - OPTIONAL. Base directory of solar flux data. Used to plot calibration bias from sun monitoring.
   * - iconpath
     - STRING
     - OPTIONAL. Base directory of ICON (NWP model) data files.



.. admonition:: Note on S3 and MySQL writing

   Every Pyrad product output (PNG or CSV) can be written to an S3 bucket,
   provided that **s3EndpointWrite** and **s3BucketWrite** are specified in
   the main configuration file and the environment variables
   **S3_KEY_WRITE** and **S3_SECRET_WRITE** contain the access key ID and
   secret access key at runtime. The Python library **boto3** must be
   installed.

   Every CSV product output can also be written to a MySQL database,
   provided that **MySQLhost**, **MySQLdatabase**, and **MySQLuser** are
   specified in the main configuration file and the environment variable
   **MYSQL_PASSWORD** contains the database password at runtime. The
   Python libraries **sqlalchemy**, **cryptography**, and **pymysql**
   must be installed.

   For each product, **s3write** and **MySQLwrite** control whether its
   output is written to S3 or MySQL. By default, **s3write** is ``True``
   and **MySQLwrite** is ``False``. Writing requires the corresponding
   connection settings and credentials.

   The MySQL table name is derived from the configured product name.
   If a product generates multiple CSV files, an index is appended to
   the table name, for example ``GC_monitoring_Hpol_1``.

   When a datetime column is detected, the first detected datetime
   column is used as the time key. Duplicate timestamps within an input
   CSV are reduced to the last occurrence. By default, rows whose
   timestamps already exist in the database are skipped.

   The optional configuration keyword **MySQLoverwrite** defaults to
   ``False``. When set to ``True`` in the product configuration file,
   rows with matching timestamps are updated and rows with new
   timestamps are inserted. All non-time columns are updated,
   including replacing existing values with SQL ``NULL`` when the
   corresponding CSV values are missing. Database rows whose
   timestamps are absent from the CSV remain unchanged. Missing or
   invalid time-key values cause the import to fail when overwriting
   is enabled.

   If no datetime column is detected, exact duplicate rows are skipped
   by comparing all columns. In this case, **MySQLoverwrite** has no
   effect.

   The optional configuration keyword **MySQLdrop** defaults to
   ``False``. When set to ``True`` in the product configuration file,
   an existing table is dropped and recreated if its column names or
   column order differ from those of the CSV. This deletes all existing
   rows in that table. Column data types are not compared by this check.
   When **MySQLdrop** is ``False``, a column mismatch prevents the import
   and raises an error.
