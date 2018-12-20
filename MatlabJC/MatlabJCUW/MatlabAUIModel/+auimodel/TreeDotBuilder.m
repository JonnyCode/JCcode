classdef TreeDotBuilder < handle
  properties
    dotString = '';
    visitor;
  end

  methods
    function buildDot(self, tree)
      import auimodel.TreeDotVisitor;

      self.visitor = TreeDotVisitor();
      tree.accept(self.visitor);
      self.visitor.finalize();
    end

    function displayDiagram(self)
      dotPath = self.tempfile('dot');

      fid = fopen(dotPath, 'w');

      iter = self.visitor.dotLines.iterator;
      while iter.hasNext
        line = iter.nextValue;
        fprintf(fid, [line '\n']);
      end

      status = fclose(fid);

      if ~ status == 0
        message = ferror(fid);
        MException('AUIModel:SystemFailed', ...
          ['Failed to write "' dotPath '": ' message]).throw();
      end

      pngPath = self.tempfile('png');
      self.runAndCheck(['/usr/local/bin/dot -Tpng -o ' pngPath ' ' dotPath]);
      self.runAndCheck(['open ' pngPath]);
    end

    function runAndCheck(self, systemCommand)
      [status, output] = system(systemCommand);

      if ~ status == 0
        MException('AUIModel:SystemFailed', ...
          ['System command "' systemCommand '" failed: ' output]).throw();
      end
    end

    function path = tempfile(self, extension)
      % HACK: ugh, use tempname instead (and check reference more often)
      path = ['/tmp/auimodel' num2str(round(1000*rand(1))) '.' extension];
    end
  end
end

