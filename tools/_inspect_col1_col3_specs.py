r"""column1/40.hsc と column3/150.hsc の Specifications を列挙。Comp Recovery が Active か確認。"""
import os, sys
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))
from units.vle.hysys.registry import HsysRegistry
from units.vle.hysys.session import HysysSession

def dump_specs(label, hsc_path):
    print(f"\n===== {label}: {hsc_path.name} =====")
    with HysysSession(hsc_path, visible=False) as sess:
        col = sess.flowsheet.Operations.Item("T-100")
        specs = col.ColumnFlowsheet.Specifications
        for i in range(specs.Count):
            s = specs.Item(i)
            try:
                name = s.Name
                ia = s.IsActive
                cv = s.CurrentValue
                gv = s.GoalValue
                cv_str = f"{cv:.4f}" if abs(cv - (-32767)) > 1 else "EMPTY"
                gv_str = f"{gv:.4f}" if abs(gv - (-32767)) > 1 else "EMPTY"
                print(f"  [{i}] {name}: Active={ia}, Cur={cv_str}, Goal={gv_str}")
            except Exception as e:
                print(f"  [{i}] error: {e}")

if __name__ == "__main__":
    reg = HsysRegistry()
    dump_specs("column1", reg.get_path("column1", 40))
    dump_specs("column3", reg.get_path("column3", 150))
