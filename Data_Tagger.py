import itertools
class Scheme():
    
    # Cryptic errors when incorrrect slash
    # Cryptic error when the number of attr given exceeds what exists in the filepath
    
    def __init__(self,  attr_list, split_list=[], ignore_attr='_',
                 val_delim='_', split_delim='-', incl_extension=False, slash='/'):
        
        
        self.attr_list = attr_list
                    
        self.split_list = split_list
        self.ignore_attr = ignore_attr 
        self.val_delim = val_delim
        self.split_delim = split_delim
        self.incl_extension = incl_extension
        self.slash = slash

        self.attrs = [attr for attr in self.attr_list if attr not in self.ignore_attr + self.slash]

        

    def __str__(self):
        return self.scheme_name

    
    
    def _pair_attr_vals(self, filename):
        
        # replace file extension dot with delim if incl_extension == True
        if self.incl_extension == True:          
            fname = filename[::-1].replace('.',self.scheme_obj.val_delim)[::-1]
        else: 
            fname = filename.rpartition('.')[0]
        
#         print(self.name)
        
        # Break the scheme list and file str into blocks based on directory
        schem_dir_blocks = [list(group) for buul, group in itertools.groupby(self.attr_list, lambda x: x == self.slash) if not buul]
        n_schem_dir_blocks = len(schem_dir_blocks)
        
        file_dir_blocks = fname.rsplit(self.slash, maxsplit=n_schem_dir_blocks)[-n_schem_dir_blocks:]
        
#         print('file_dir_blocks')
#         print(file_dir_blocks)
    
        # Break scheme dir blocks further up by the '*'
        schem_star_blocks = []
        for block in schem_dir_blocks:
            if block.count('*') > 1:
                raise Exception('More than 1 * cannot be used per directory string')
            
            if block[0] == '*':
                block.remove('*')
                schem_star_blocks.append([[],block])
            elif block[-1] == '*':
                block.remove('*')
                schem_star_blocks.append([block,[]])
            else:
                sublock = [list(group) for buul, group in itertools.groupby(block, lambda x: x == "*") if not buul]
                schem_star_blocks.append(sublock)
        
#         print('scheme_star_blocks:')
#         print(schem_star_blocks)
        
        # assign val-attr
        attr_val_list = []
        for i in range(0, n_schem_dir_blocks):
            
            file_dir_val_list = file_dir_blocks[i].split(self.val_delim)
            n_star_blocks = len(schem_star_blocks[i])
            
            # If there are no stars in this dir (and so only one block):
            if n_star_blocks == 1:
                self._check_lengths(fname, schem_star_blocks[i][0], file_dir_val_list)
                va = list(zip(schem_star_blocks[i][0], file_dir_val_list))
                attr_val_list.append(va)
            
            # If there is a star (as so two blocks that were divided by that star)
            elif n_star_blocks == 2:
                # Slices. The first scheme star block zips from the front, the second from the back
                # The vals from the file string are sliced by the number of attr in the scheme star blocks
                
                if (len(schem_star_blocks[i][0]) + len(schem_star_blocks[i][1])) > len(file_dir_val_list):
                    print(file_dir_val_list)
                    print(schem_star_blocks[i][0])
                    print(schem_star_blocks[i][1])
                    raise Exception('ERROR: There are more attributes than values!!!')
                
                
                va = list(zip(schem_star_blocks[i][0], file_dir_val_list[:len(schem_star_blocks[i][0])]))
                attr_val_list.append(va)
                
                va = list(zip(schem_star_blocks[i][1], file_dir_val_list[-(len(schem_star_blocks[i][1])):]))
                attr_val_list.append(va)
                
            else:
                raise Exception('ERROR: More than one star has been used in a directory!!!')
    
        # flatten list and remove ignored attr and split attr
        flatl = []     
        for sublist in attr_val_list:
            for av in sublist:
                av = list(av)
#                 print(av)
                if av[0] == self.ignore_attr:
                    continue
                if av[0] in self.split_list:
                    av[1] = av[1].split(self.split_delim)[-1]
                flatl.append(av)
                
        # convert flatl to dict
        attr_val_dict = {i[0]:i[1] for i in flatl}
        
        return attr_val_dict
        
    def _check_lengths(self, file, attrs, vals):
        if len(vals) != len(attrs):
            print(f'filename: {file}')
            print('The attributes:')
            print(attrs)
            print('are being paired to the values:')
            print(vals)
            print('')
            raise Exception('ERROR: Number of values does not match number of attributes!!!')


class Base_Manager():

    def __init__(self, ):


        self.objs = []
        self.attr = []
        self.file_attr = {}



    # BUGS:
    #

    # Potential issues:
    # What if the order of the attr between file names is different?
    # What if its a first upload and you want several obj to one file?
    # No common attr also triggers the non unique file exception - which it doesn't excplicitly relate to
    # The same error message could probably be less ambiguous as well
    # If a unique attr is  ignored in the scheme, two files with the same attr maybe present and cause ambiguity
    # Multiple files to an obj might be useful in some situations

    # ToDo:
    # If uploaded onto same attrname, perform some sort of merge operation
    # Aggregate upload of data, or aggregate access in some way


    def add_analysis(self, file_list, scheme=None, uploader=None,):

        """
        This function works out how newly uploaded files should be paired with existing objects.

        It first finds what attributes the objects and the new files have in common.

        The data, and new attr, will be added to the objects which share the same attribute values.

        The intersection of attributes between the objs and the new files must only have a single file associated
        with it. This is so that there is no conflict in what the data is assigned to an attr
        """

        # add checks to make sure things aren't none and also save

        if len(file_list) == 0:
            raise Exception("The file list was empty")
        if type(file_list) != list:
            raise Exception('file_list must be a list of strs')

        # check that the attribute from the uploader has not been used before, and no data will be overwritten
        if uploader.attrname in self.attr:
            raise Exception(f'The attribute {uploader.attrname} has already been used by this Manager obj')


        # Find on which objs the file contents should be stored
        # And also add new attr to those obj
        if len(self.objs) == 0:
            file_obj_pairings = self._init_analysis(file_list, scheme, uploader)
        else:
            file_obj_pairings = self._pair_objs_with_files(file_list, scheme, )

#         print(file_obj_pairings)
        # Run the uploader - Put the file contents into the objects
        for file, objs in file_obj_pairings.items():
            uploader._upload(file, objs)

        # Make the colection of the new data accesible from the manager object
        setattr(self, uploader.attrname, self.get_values(uploader.attrname) )


    def _init_analysis(self, file_list, scheme, uploader):


        """
        Uploads files and creates object representations of them. This function is called
        when there are no existing objects already that the new files must be linked with.
        """

        file_obj_pairings = {}

        for file in file_list:

            obj = Base_Obj(file, scheme)
            self.file_attr[file] = obj.attr_vals

            for attr in obj._get_attr():
                if attr not in self.attr:
                    self.attr.append(attr)

            file_obj_pairings[file] = [obj]
            self.objs.append(obj)

        return file_obj_pairings




    def _pair_objs_with_files(self, file_list, scheme,):

        """
        Determines which pre existing objects correspond to the newly specified files on the basis of
        which attributes the objects and files have in common.

        Attributes unique to the new files are also added to the corresponding objects after pairing.

        """


        ##########################################
        ## Pair objects with files
        # Change _pari_attr_val to return a dictionary
        # dictionary of dictionaries
        attr_val_dicts = {file: scheme._pair_attr_vals(file) for file in file_list}
        self.file_attr.update(attr_val_dicts)

        # First see what attr are in common between the existing obj and the files to be added
#         common_attr = self.attr.intersection(set(scheme.attr_list))
#         common_attr = self.attr.keys() & scheme.attr_list

        common_attr = [attr for attr in scheme.attr_list if attr in self.attr ]


        # Check that the values of the common attr are enough to define each file as unique
        # If the common values already exist in the set, they will not be added, and
        # there will be a different number of files and common value tuples in the set.

        # also gen a dicitonary of common values: filename, to help pair objects later
        common_attr_sets = set()
        common_value_filenames = {}

        for filename, av_dict in attr_val_dicts.items():

            common_values = tuple([av_dict[attr] for attr in av_dict if attr in common_attr])
            common_attr_sets.add(common_values)

            common_value_filenames[common_values] = filename

#         print('com val filenames')
#         print(common_value_filenames)

        if len(common_attr_sets) != len(file_list):
            print()
            print(common_attr_sets)
            print()
            print(file_list)

            raise Exception('The attributes given to the file are not enough to uniquely define each file')


        # break the objects by the common attr.
        # The file with the corresponding values of those attributes is paired with the objects
        # need to associate the value tuple with the filename
        # sort

        get_common_values = lambda obj: tuple([getattr(obj, attr) for attr in common_attr])

        self.objs.sort(key=get_common_values)

#         print([get_common_values(o) for o in self.objs])


        obj_file_pairings = {}


        for common_values, objs_gen in itertools.groupby(self.objs, key=get_common_values):

            objs = list(objs_gen)
            filename = common_value_filenames[common_values]
            obj_file_pairings[filename] = objs



            # update obj with all attr from the new file
            for obj in objs:
                obj._label(**self.file_attr[filename])

        # update attr list
        self._update_attr()


        return obj_file_pairings




    def _update_attr(self,):
        """ Assumes that all objects have the same attr """
        self.attr = self.objs[0]._get_attr()




    def __iter__(self, ):
        return (o for o in self.objs)

    def get_values(self, *attrs, unique=False):

        """
        Returns all values of an inputted attribute across the all objects
        Unique (Bool): Return only unique values


        """

        out = []
        for attr in attrs:
            vals = []
            for obj in self.objs:
                    if hasattr(obj, attr):
                        val = getattr(obj, attr)
                        vals.append(val)
                    else:
                        vals.append(np.nan)

            out.append(vals)


        if len(attrs) == 1:
            out = out[0]

        else:
            out = list(zip(*[vals for vals in out]))

        if unique == True:
            out = list(set(out))

        return out


    def get_objects(self, attr_val):

        """
        Returns all objects that match the given attribute - value pairs, specified as a dict.
        Values can be a list if multiple values of an attribute are acceptable

        """

        # If attr_val is empty, all objects are returned - which I don't think is hte behaviour that i want 

        attr_val_to_search = attr_val.copy()

        out = []

        for attr, vals in attr_val_to_search.items():
            if not isinstance(vals, list):
                attr_val_to_search[attr] = [vals]

        for obj in self.objs:
            attr_tf = []
            for attr, vals in attr_val_to_search.items():

                if any([getattr(obj, attr) == val for val in vals]):
                    attr_tf.append(True)
                else:
                    attr_tf.append(False)

            if all(attr_tf):
                out.append(obj)


        return out



class Base_Obj():

    def __init__(self, fname, scheme ):
        self.filename = fname
        self.scheme = scheme

        self._added_attr = []

        self.attr_vals = self.scheme._pair_attr_vals(self.filename)

        self._label(**self.attr_vals)



    def _label(self, **kwargs):


        for attr, val in kwargs.items():
            setattr(self, attr, val)
            self.attr_vals.update({attr:val})



    def _get_attr(self, ):
        return list(self.attr_vals.keys())


    def __repr__(self,):
        return "Data Object < " + ", ".join([ f'{attr}: {val}' for attr, val in self.attr_vals.items()])+" >"


class Uploader():


    def __init__(self, attrname, spliton=None, splitfunc=None, uploadfunc=None, **kwargs):

        """
        Uploader class stores functions and other kwargs for the purposes of uploading data from a
        file, and assigning it to objects in the manner desired.

        kwargs:

        attrname (str): Name given to the attribute that stores the uploaded file's data on an object

        spliton (str): Several objects may be associated with the information in a file. If None,
                        each object is given all the data from the uploaded file. If it is a str, it
                        must correspond to an existing attribute on the objects. The value of the
                        attribute of that object will be passed to kwarg splitfunc, which will split the data
                        based on said value and store the corresponding data on the attribute attrname of
                        that object.

        splitfunc (func): function object which should only take the positional arguments, data_structure and
                        getattr(obj, spliton)

        uploadfunc (func): Function which uploads the file into some data structure. The first positional
                        argument must be a filepath str and then kwargs.

        kwargs (dict): Kwargs to be passed to uploadfunc.
        """


        self.attrname = attrname
        self.spliton = spliton
        self.splitfunc = splitfunc
        self.uploadfunc = uploadfunc
        self.kwargs = kwargs

        if (spliton == None) != (splitfunc == None):
            raise Exception("Either spliton or splitfunc kwarg has been defined and not the other")



    def _upload(self, file, objs, ):

        data = self.uploadfunc(file, **self.kwargs)
#         print(data)
        for obj in objs:

            if self.spliton == None:
                try:
                    setattr(obj, self.attrname, copy(data))
                except:
                    setattr(obj, self.attrname, data)
            else:
                data_slice = self.splitfunc(data, getattr(obj, self.spliton))
                setattr(obj, self.attrname, data_slice)



