import os
import logging

def compile(args):

    logger = logging.getLogger("USCM")

    logger.info("Compiling UAMMD launcher ...")

    try:

        import subprocess
        from setuptools_cuda_cpp import find_cuda_home_path

        try:

            UAMMD_PATH = os.environ['UAMMD_PATH']

        except:

            logger.error("Environment variables UAMMD_STRUCTURED_PATH")
            raise RuntimeError("Environment variables not found")

        try:

            UAMMD_STRUCTURED_PATH = os.environ['UAMMD_STRUCTURED_PATH']

        except:

            logger.error("Environment variables UAMMD_PATH")
            raise RuntimeError("Environment variables not found")

        INCLUDERS = ['-I'+UAMMD_PATH+'/src/',
                     '-I'+UAMMD_PATH+'/src/third_party/',
                     '-I'+UAMMD_STRUCTURED_PATH+'/']

        try:
            CUDA_PATH = find_cuda_home_path()
        except:
            raise RuntimeError('Can not find CUDA_HOME path')

        if CUDA_PATH is None:
            raise RuntimeError('Can not find CUDA_HOME path')

        INCLUDERS.append('-I'+str(CUDA_PATH)+'/include/')

        INCLUDERS = ' '.join(INCLUDERS)

        FLAGS = ' '.join(['--expt-relaxed-constexpr',
                          '--expt-extended-lambda',
                          '-std=c++14',
                          '-O3',
                          '-DUAMMD_EXTENSIONS',
                          '-DMAXLOGLEVEL=5',
                          #'-Xcompiler=\"-O3 -march=native -fPIC\"',
                          '-Xcompiler=\"-O3 -fPIC\"',
                          '-ccbin=g++',
                          '-w'])

        LIBRARIES = ' '.join(['-lcufft',
                              '-llapacke',
                              '-lcublas',
                              '-lblas',
                              '-lcurand',
                              '-lcusolver',
                              '-lcusparse',
                              '-lstdc++fs'])
        try:
            CUDA = os.path.join(find_cuda_home_path(),'bin','nvcc')
        except:
            raise RuntimeError('Can not find CUDA_HOME path')

        if args.ccbin:
            ccbin = args.ccbin[0]
            FLAGS = [ f for f in FLAGS.split(" ") if not f.startswith("-ccbin") ]
            FLAGS.append("-ccbin="+ccbin)
            FLAGS = " ".join(FLAGS)

        if args.mkl:
            FLAGS += " -DUSE_MKL"
            #Remove lapacke and blas from libraries
            LIBRARIES = [ l for l in LIBRARIES.split(" ") if not l.startswith("-llapacke") and not l.startswith("-lblas") ]
            LIBRARIES += "-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core".split(" ")
            LIBRARIES += "-L/opt/intel/lib/intel64 -liomp5 -lpthread -lm -ldl".split(" ")
            LIBRARIES = " ".join(LIBRARIES)

        if not args.arch:

            try:
                NVCC_PATH = os.path.join(find_cuda_home_path(), 'bin', 'nvcc')

                gencode = []

                # Get the current list gpu capturing nvcc output
                p = subprocess.Popen([NVCC_PATH, '--list-gpu-code'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()

                archs = out.decode('utf-8').split('\n')[0:-1]
                archs = [arch.split("_")[-1] for arch in archs]

                for arch in archs:
                    gencode += ['-gencode',f'arch=compute_{arch},code=sm_{arch}']

                ARCH = ' '.join(gencode)

            except:
                raise RuntimeError('Could not run nvcc')

        else:
            gencode = []
            for a in args.arch:
                gencode += ['-gencode',f'arch=compute_{a},code=sm_{a}']
            ARCH = ' '.join(gencode)

        if args.verbose:
            FLAGS=FLAGS.replace("DMAXLOGLEVEL=5","DMAXLOGLEVEL=7")
            FLAGS=FLAGS.replace("-w","")

        if args.double:
            FLAGS += " -DDOUBLE_PRECISION"

        #####################################

        UAMMD_PATH = os.environ["UAMMD_PATH"]

        if args.folder:
            folder = args.folder[0]
        else:
            folder = UAMMD_PATH+'/bin/'

        if not os.path.exists(folder):
            self.logger.error("Folder {} does not exist".format(folder))
            raise RuntimeError("Folder does not exist")

        if args.name:
            name = args.name[0]
        else:
            name = 'UAMMDlauncher'

        OUTPUT  = os.path.join(folder,name)
        COMPILE = " ".join([CUDA,FLAGS,ARCH,UAMMD_PATH+"/launcher/UAMMDlauncher.cu",INCLUDERS,LIBRARIES,"-o",OUTPUT])

        logger.info(COMPILE)

        subprocess.run([COMPILE],shell=True)

    except Exception as e:
        logger.error("Error compiling UAMMD launcher: "+str(e))
        raise Exception("Compilation error")
