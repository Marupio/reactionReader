    Foam::Info<< "Create admTime\n" << Foam::endl;

    Foam::admTime runTime
    (
        Foam::admTime::controlDictName,
        args.rootPath(),
        args.caseName()
    );
